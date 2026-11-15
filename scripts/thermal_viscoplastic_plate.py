"""Quarter plate: prescribed thermal ramp and hold with optional elastic softening."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.nonlinear_solver import NonlinearModel, solve_increment
from immersed_elasticity.solver import structured_mesh
from immersed_elasticity.boundary_conditions import thermal_constraints
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.viscoplasticity import J2Material
from immersed_elasticity.thermal_settings import thermal_settings
from immersed_elasticity.thermal_plots import save_thermal_plots

def main(argv=None, *, initial_settings=None, output=None, overwrite=False):
    parser = argparse.ArgumentParser(description=__doc__)
    try:
        initial = thermal_settings(initial_settings)
    except (ValueError, TypeError) as error:
        parser.error(str(error))
    parser.add_argument("--mesh", type=int, default=None)
    parser.add_argument("--depth", type=int, default=initial["quadrature"]["depth"])
    parser.add_argument("--alpha", type=float, default=initial["quadrature"]["alpha"])
    parser.add_argument("--method", choices=("quadtree", "moment-fitting"), default=initial["quadrature"]["method"],
                        help="Mechanical integration only; heat always uses quadtree.")
    parser.add_argument("--steps-per-second", type=int, default=initial["loading"]["steps_per_second"])
    parser.add_argument("--right-temperature", type=float, default=initial["loading"]["right_temperature"])
    parser.add_argument("--elasticity", choices=("constant", "paper-reference", "paper-temperature"),
                        default=initial["properties"]["elasticity"], help="Use original constants, paper fits frozen at 20 C, or local-temperature fits.")
    parser.add_argument("--hardening", choices=("constant", "paper-reference", "paper-temperature"),
                        default=initial["properties"]["hardening"], help="Use original Q/b, paper Q/b at 20 C, or temperature-dependent Q.")
    parser.add_argument("--supports", choices=("symmetry", "minimal"), default=initial["supports"])
    parser.add_argument("--output", type=Path, default=Path(output) if output is not None else Path("results/walkthrough/v4-thermal-viscoplastic-plate"))
    parser.add_argument("--overwrite", action="store_true", default=overwrite)
    args = parser.parse_args(argv)
    if args.mesh is not None:
        initial["geometry"].update(nx=args.mesh, ny=args.mesh)
    initial["quadrature"].update(depth=args.depth, alpha=args.alpha, method=args.method)
    initial["loading"].update(steps_per_second=args.steps_per_second, right_temperature=args.right_temperature)
    initial["properties"].update(elasticity=args.elasticity, hardening=args.hardening)
    initial["supports"] = args.supports
    try:
        inputs = thermal_settings(initial)
    except (ValueError, TypeError) as error:
        parser.error(str(error))
    geometry, loading = inputs["geometry"], inputs["loading"]
    material = J2Material(**inputs["material"])
    thickness, conductivity = inputs["thickness"], inputs["conductivity"]
    reference_temperature = loading["reference_temperature"]
    expansion_coefficient = inputs["expansion_coefficient"]
    heating_time, hold_time = loading["heating_time"], loading["hold_time"]
    if args.output.exists() and (not args.output.is_dir() or
                                 (any(args.output.iterdir()) and not args.overwrite)):
        parser.error("Choose an empty output directory or use --overwrite.")
    mesh = structured_mesh(**geometry)
    x, y = mesh.nodes.T
    left = np.flatnonzero(np.isclose(x, 0.) & (y >= mesh.radius))
    right = np.flatnonzero(np.isclose(x, mesh.length))
    temperatures = {int(i): loading["left_temperature"] for i in left}
    temperatures.update({int(i): args.right_temperature for i in right})
    heat = solve_stationary_heat(mesh, conductivity, thickness, temperatures,
                                depth=args.depth, alpha=args.alpha)
    model = NonlinearModel(mesh, material, thickness, depth=args.depth, alpha=args.alpha,
                           elasticity=args.elasticity, hardening=args.hardening, method=args.method)
    constraints = thermal_constraints(mesh, args.supports)
    state = model.initial_state()
    history = [[0., reference_temperature, 0., 0., 0., 0., 0]]
    heating_steps = int(np.ceil(heating_time*args.steps_per_second))
    hold_steps = int(np.ceil(hold_time*args.steps_per_second))
    times = np.r_[0., np.arange(1, heating_steps+1)*(heating_time/heating_steps),
                  heating_time+np.arange(1, hold_steps+1)*(hold_time/hold_steps) if hold_steps else []]
    maximum_residual = maximum_imbalance = 0.
    maximum_iterations = 0
    print(f"Active cells: {len(mesh.cells)}; DOFs: {model.ndofs}; integration points: {model.point_count}", flush=True)
    print(f"Mechanical method: {args.method}; physical history points: {np.count_nonzero(model.physical)}; thermal method: quadtree", flush=True)
    print(f"0-{heating_time:g} s: right edge {reference_temperature:g} -> {args.right_temperature:g} C; "
          f"hold to {heating_time+hold_time:g} s; maximum dt={np.max(np.diff(times)):g} s", flush=True)
    print(f"Left edge ramps from {reference_temperature:g} to {loading['left_temperature']:g} C; supports: {args.supports}.", flush=True)
    print(f"Elasticity: {args.elasticity}; hardening: {args.hardening}; initial yield, viscosity, and expansion coefficient fixed.", flush=True)
    start = perf_counter()
    progress_interval = max(1, (len(times)-1)//10)
    for step in range(1, len(times)):
        time, dt = times[step], times[step]-times[step-1]
        scale = min(time/heating_time, 1.)
        temperature = reference_temperature + scale*(heat.temperature-reference_temperature)
        result = solve_increment(model, state, constraints, dt, nodal_temperature=temperature,
                                 reference_temperature=reference_temperature,
                                 expansion_coefficient=expansion_coefficient)
        state = result.state  # History is committed only after global convergence.
        stress = result.stress[model.physical]
        von_mises = np.sqrt(stress[:, 0]**2-stress[:, 0]*stress[:, 1]+stress[:, 1]**2+3*stress[:, 2]**2)
        plastic = np.array([h.accumulated_strain for h in state.history])[model.physical]
        ux = state.displacement.reshape(-1, 2)[right, 0].mean()
        imbalance = np.max(np.abs(result.residual.reshape(-1, 2).sum(axis=0)))
        maximum_residual = max(maximum_residual, result.free_residual_norm)
        maximum_imbalance = max(maximum_imbalance, imbalance)
        maximum_iterations = max(maximum_iterations, result.iterations)
        history.append([time, reference_temperature+scale*(args.right_temperature-reference_temperature),
                        ux, von_mises.max(), plastic.max(), result.free_residual_norm, result.iterations])
        if step in (heating_steps, len(times)-1):
            print(f"t={time:g} s: mean right ux={ux:.9f} mm; max von Mises={von_mises.max():.6f} MPa; "
                  f"max accumulated strain={plastic.max():.8f}", flush=True)
        elif step % progress_interval == 0:
            print(f"Progress: {step}/{len(times)-1} increments; t={time:g} s", flush=True)
    print(f"Maximum free residual: {maximum_residual:.3e} N; force imbalance: {maximum_imbalance:.3e} N")
    solve_seconds = perf_counter()-start
    print(f"Maximum Newton corrections: {maximum_iterations}; solve time: {solve_seconds:.3f} s")
    print(f"Plane-stress updates: {model.plane_stress_updates}; global evaluations: {model.evaluations}")
    print("Prescribed equilibrium temperature snapshots; no plastic heating.")
    point_temperature = np.concatenate([shape @ temperature[cell] for cell, shape in model.temperature_interpolation])[model.physical]
    peak = int(np.argmax(plastic))
    location = model.points[model.physical][peak]
    print(f"Maximum-plastic-strain sample: x={location[0]:.6f}, y={location[1]:.6f} mm; T={point_temperature[peak]:.3f} C")
    args.output.mkdir(parents=True, exist_ok=True)
    settings = dict(geometry=geometry,
                    material=asdict(material), elasticity=args.elasticity, hardening=args.hardening,
                    reference_elastic_material=asdict(model.material),
                    conductivity_W_per_mm_K=conductivity, thickness_mm=thickness,
                    reference_temperature_C=reference_temperature, left_temperature_C=loading["left_temperature"],
                    right_temperature_C=args.right_temperature, expansion_coefficient_per_K=expansion_coefficient,
                    depth=args.depth, alpha=args.alpha, steps_per_second=args.steps_per_second,
                    mechanical_method=args.method, thermal_method="quadtree",
                    heating_time_s=heating_time, hold_time_s=hold_time, supports=args.supports, inputs=inputs)
    (args.output / "settings.json").write_text(json.dumps(settings, indent=2)+"\n", encoding="utf-8")
    metrics = dict(integration_points=model.point_count, physical_history_points=int(np.count_nonzero(model.physical)),
                   plane_stress_updates=model.plane_stress_updates, global_evaluations=model.evaluations,
                   maximum_newton_corrections=maximum_iterations, solve_seconds=solve_seconds,
                   maximum_free_residual_N=maximum_residual, maximum_force_imbalance_N=maximum_imbalance,
                   fallback_cells=model.fallback_cells)
    (args.output / "metrics.json").write_text(json.dumps(metrics, indent=2)+"\n", encoding="utf-8")
    history = np.asarray(history)
    np.savetxt(args.output / "history.csv", history, delimiter=",", comments="",
               header="time_s,right_temperature_C,mean_right_ux_mm,max_von_mises_MPa,max_accumulated_strain,free_residual_N,newton_corrections")
    np.savetxt(args.output / "nodes.csv", np.column_stack((mesh.nodes, temperature, state.displacement.reshape(-1, 2))),
               delimiter=",", comments="", header="x_mm,y_mm,temperature_C,ux_mm,uy_mm")
    np.savetxt(args.output / "physical_samples.csv", np.column_stack((model.points[model.physical], stress, von_mises, plastic, point_temperature)),
               delimiter=",", comments="", header="x_mm,y_mm,sigma_xx_MPa,sigma_yy_MPa,tau_xy_MPa,von_mises_MPa,accumulated_strain,temperature_C")
    save_thermal_plots(args.output, history, model.points[model.physical], point_temperature,
                       von_mises, plastic, mesh, loading)
    summary = (
        "CutCell Mechanics - run summary\n\n"
        f"Mechanical supports: {args.supports}\n"
        "symmetry = left x rollers and bottom y rollers; minimal = rigid motion removed only.\n"
        f"Mesh: {geometry['nx']} x {geometry['ny']}; method: {args.method}\n"
        f"Elasticity: {args.elasticity}; hardening: {args.hardening}\n"
        f"Final mean right-edge x displacement: {ux:.9f} mm\n"
        f"Final sampled peak von Mises stress: {von_mises.max():.6f} MPa\n"
        f"Final maximum accumulated plastic strain: {plastic.max():.8g}\n"
        f"Temperature at that plastic-strain sample: {point_temperature[peak]:.3f} C\n"
        f"Maximum free equilibrium residual: {maximum_residual:.3e} N\n"
        f"Solve time: {solve_seconds:.3f} s (single run, not a performance benchmark)\n\n"
        "Open history.png to see the heating schedule and mechanical response.\n"
        "The shaded interval is the temperature hold. Plastic strain records permanent\n"
        "deformation history; zero means no plastic accumulation in this run.\n"
        "Open fields.png to locate the final temperature, stress, and plastic strain.\n"
        "Dots are physical integration points, not a finer displacement mesh.\n"
        "Colour limits vary between fields/runs. Compare numerical values as well.\n\n"
        "This is a coarse demonstration unless resolution has been verified separately.\n"
        "The restrained-heating hotspot can be near the cold supported corner.\n"
        "Low residual confirms equilibrium, not mesh convergence or physical validation.\n"
        "Temperature snapshots are stationary; no transient heat or plastic heating.\n"
        "See docs/GETTING_STARTED.md and docs/MODEL.md for interpretation and limits.\n")
    (args.output / "summary.txt").write_text(summary, encoding="utf-8")
    print(f"Saved results to {args.output.resolve()}\nOpen summary.txt, history.png, and fields.png first.")
    return args.output


if __name__ == "__main__":
    main()
