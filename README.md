# Claude Scientific Skills

Claude Scientific Skills is a local collection of 141 scientific skills for people who want Claude to use domain-specific workflows, tools, and references without manually assembling prompts for each task.

## Documentation

- [docs/SKILLS_INDEX.md](docs/SKILLS_INDEX.md): Alphabetical list of skills with one-line descriptions and direct links to each `SKILL.md` file.
- [docs/scientific-skills.md](docs/scientific-skills.md): Detailed catalog of skill coverage by category.
- [docs/examples.md](docs/examples.md): Example prompts and usage patterns.
- [docs/open-source-sponsors.md](docs/open-source-sponsors.md): Upstream projects that power many skills.

## Quick start

Primary path: discover and open skills directly from this repo.

```bash
# List all skill folders that contain SKILL.md
find scientific-skills -mindepth 2 -maxdepth 2 -name SKILL.md | sort

# Read one skill
sed -n '1,200p' scientific-skills/rdkit/SKILL.md
```

If you use Claude Code plugin workflows, install and manage this repository from your client and use the skills index above to pick relevant skills faster.

## License

- [LICENSE.md](LICENSE.md)
