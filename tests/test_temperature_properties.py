import numpy as np
import pytest

from immersed_elasticity.temperature_properties import with_temperature_dependent_elasticity
from immersed_elasticity.viscoplasticity import J2Material
from immersed_elasticity.temperature_properties import with_temperature_dependent_hardening


def test_hardening_matches_paper_reference_curve_and_preserves_other_properties():
    reference = J2Material(initial_yield=302.3605, relaxation_time=100)
    material = with_temperature_dependent_hardening(reference, 20)
    np.testing.assert_allclose(material.saturation_increment, 206.3425, rtol=1e-12)
    for p in (0., .02, .1):
        expected = 302.3605 + 206.3425 * (1 - np.exp(-p / .02784959))
        np.testing.assert_allclose(material.yield_stress(p), expected, rtol=1e-12)
    for name in ("young", "poisson", "initial_yield", "relaxation_time", "rate_exponent"):
        assert getattr(material, name) == getattr(reference, name)


def test_hardening_softens_at_example_temperatures_without_changing_initial_yield():
    materials = [with_temperature_dependent_hardening(J2Material(), t) for t in (20, 400, 700)]
    assert np.all(np.diff([m.yield_stress(.02) for m in materials]) < 0)
    assert all(m.yield_stress(0) == 250. for m in materials)


@pytest.mark.parametrize("temperature", [-1, 1501, np.nan, np.inf])
def test_hardening_rejects_invalid_temperature(temperature):
    with pytest.raises(ValueError, match="0-1500"):
        with_temperature_dependent_hardening(J2Material(), temperature)


def test_reference_temperature_recovers_paper_bulk_and_shear_moduli():
    reference = J2Material(relaxation_time=100)
    material = with_temperature_dependent_elasticity(reference, 20)
    np.testing.assert_allclose(material.bulk_modulus, 163195.2, rtol=1e-12)
    np.testing.assert_allclose(material.shear_modulus, 80891.53, rtol=1e-12)
    assert material.initial_yield == reference.initial_yield
    assert material.saturation_increment == reference.saturation_increment
    assert material.hardening_rate == reference.hardening_rate
    assert material.relaxation_time == reference.relaxation_time
    assert material.rate_exponent == reference.rate_exponent
    assert reference.young == 210000.


def test_elastic_stiffness_decreases_at_the_example_temperatures():
    materials = [with_temperature_dependent_elasticity(J2Material(), temperature)
                 for temperature in (20, 400, 700)]
    assert np.all(np.diff([material.young for material in materials]) < 0)
    assert np.all(np.diff([material.shear_modulus for material in materials]) < 0)
    assert all(material.bulk_modulus > 0 and 0 < material.poisson < .5 for material in materials)


@pytest.mark.parametrize("temperature", [-1, 1501, np.nan, np.inf])
def test_fit_rejects_out_of_range_or_nonfinite_temperature(temperature):
    with pytest.raises(ValueError, match="0-1500"):
        with_temperature_dependent_elasticity(J2Material(), temperature)
