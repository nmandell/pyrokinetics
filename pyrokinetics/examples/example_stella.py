import pyrokinetics

stella_template = pyrokinetics.template_dir / "input.stella"

#print("This is the stella_template: ", stella_template)

pyro = pyrokinetics.Pyro(gk_file=stella_template, gk_code="stella")

flags = {
    "stella_diagnostics_knobs": {
        "write_phi_vs_time": True,
        "write_omega": True,
    },
}

pyro.add_flags(flags)
pyro.write_gk_file(file_name="step.stella")

pyro.gk_code = "GS2"
pyro.write_gk_file(file_name="step.gs2")

pyro.gk_code = "CGYRO"
pyro.write_gk_file(file_name="step.cgyro")

pyro.gk_code = "GENE"
pyro.write_gk_file(file_name="step.gene")
