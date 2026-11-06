# Guidance for AI agents

abaqus_ufl is a formulation framework: it generates Abaqus UMAT and UEL Fortran
from declared theory.

Start here:

1. [`ai_skills/abaqus-ufl/SKILL.md`](ai_skills/abaqus-ufl/SKILL.md): which resource to use for a task, and the checks a
   claim needs. Read it before changing a model, an element or a generator.
2. [`docs/lessons/README.md`](docs/lessons/README.md): find a lesson by your task or by the symptom you see.

Repository rules:

- Work in this checkout and record its revision with any result. Other copies of abaqus_ufl or of its skill on the
  machine may be out of date.
- Change the declaration or the generator, then regenerate; do not hand-edit generated Fortran.
- Run `python -m pytest tests -q` (CI runs the same suite without Abaqus). After editing guidance files, run
  `python tools/check_guidance.py`.
- Keep private project names, paths and data out of public files.
