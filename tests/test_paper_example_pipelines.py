"""Regression checks for public paper-example generation entry points."""

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_corrosion_generation_uses_public_layout(tmp_path, monkeypatch):
    module = _load_module(
        "_test_phasefield_corrosion_cui",
        ROOT / "paper_examples" / "phasefield_corrosion"
        / "phasefield_corrosion_cui.py",
    )
    package = tmp_path / "phasefield_corrosion"
    generated = package / "generated"
    generated.mkdir(parents=True)
    monkeypatch.setattr(
        module, "__file__", str(package / "phasefield_corrosion_cui.py")
    )

    main_path, _ = module.generate()
    diag_path, _ = module.generate_diag()

    assert Path(main_path).parent == generated
    assert Path(diag_path).parent == generated
    for path in (Path(main_path), Path(diag_path)):
        text = path.read_text()
        state_write = text.index(
            "        SVARS(LOC + 15) = DBLE(ei_new_z)\n"
        )
        tangent = text.index("C       CS tangent engine", state_write)
        state_block = text[state_write:tangent]
        assert "Store fields for duplicate CPE8R" in state_block
        assert state_block.rstrip().endswith("END IF")
        assert "SUBROUTINE UVARM" in text


def test_gel_build_skips_missing_external_mesh(tmp_path, monkeypatch, capsys):
    module = _load_module(
        "_test_gel_bilayer_build",
        ROOT / "paper_examples" / "gel_bilayer" / "build.py",
    )
    missing_seed = tmp_path / "gelUEL" / "SwellInducedBending.inp"
    converted = []

    class PassingProblem:
        def verify(self, verbose=False):
            return True

    def fake_generate_uel(problem, output, **kwargs):
        Path(output).write_text("generated UEL\n")

    def fake_material_include(path):
        Path(path).write_text("generated properties\n")

    monkeypatch.setattr(module, "HERE", tmp_path)
    monkeypatch.setattr(module, "SOURCE_INP", missing_seed)
    monkeypatch.setattr(module, "ChesterAnandUPMuQuad8", PassingProblem)
    monkeypatch.setattr(module.au, "generate_uel", fake_generate_uel)
    monkeypatch.setattr(module, "_write_material_include", fake_material_include)
    monkeypatch.setattr(
        module, "_write_converted_input", lambda path: converted.append(path)
    )

    module.main()

    assert (tmp_path / "chester_anand_upmu_quad8.for").exists()
    assert (tmp_path / "ElasticGelProps_upmu_quad8.inp").exists()
    assert converted == []
    output = capsys.readouterr().out
    assert "Skipped bilayer deck regeneration" in output
    assert "abaqus/SwellInducedBending_upmu_quad8.inp" in output
