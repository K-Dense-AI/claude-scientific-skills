"""Pipeline orchestrator — chat, pipeline, and single-agent modes with Skill Registry MCP."""

from __future__ import annotations

from claude_agent_sdk import ClaudeSDKClient, ClaudeAgentOptions, query

from research_agent.agents.definitions import build_agents
from research_agent.agents.prompts import load_prompt
from research_agent.config import AgentConfig
from research_agent.skills.mcp_server import create_skill_registry_server
from research_agent.utils.message_handler import handle_message


async def run_pipeline(task_description: str, phases: list[str] | None, config: AgentConfig) -> None:
    """Run the P1-P5 pipeline by sequencing domain agent calls."""
    agents = build_agents(config)
    orchestrator_prompt = load_prompt("orchestrator.md")
    skill_server = create_skill_registry_server()

    # Build phase instruction block
    if phases:
        selected = [p for p in config.pipeline_phases if p["phase"] in phases]
    else:
        selected = config.pipeline_phases

    phase_instructions = []
    for p in selected:
        phase_instructions.append(
            f"{p['phase']} ({p['label']}): {p['domain']} 에이전트에 위임.\n"
            f"  지시: {p['instruction']}"
        )
    phase_block = "\n".join(phase_instructions)

    options = ClaudeAgentOptions(
        system_prompt=orchestrator_prompt,
        allowed_tools=["Read", "Write", "Glob", "Grep", "Bash", "Task"],
        permission_mode="acceptEdits",
        agents=agents,
        model="sonnet",
        cwd=str(config.default_cwd),
        mcp_servers={"skill-registry": skill_server},
    )

    pipeline_prompt = (
        f"다음 과제에 대해 연구 파이프라인을 실행하라:\n\n"
        f"{task_description}\n\n"
        f"실행할 단계:\n{phase_block}\n\n"
        f"각 단계에서 지정된 도메인 에이전트를 사용하고, "
        f"산출물은 파일로 저장하여 다음 단계가 참조할 수 있도록 하라."
    )

    async with ClaudeSDKClient(options=options) as client:
        await client.query(pipeline_prompt)
        async for msg in client.receive_response():
            handle_message(msg)


async def run_chat(config: AgentConfig) -> None:
    """Interactive chat mode — orchestrator routes to domain agents."""
    agents = build_agents(config)
    orchestrator_prompt = load_prompt("orchestrator.md")
    skill_server = create_skill_registry_server()

    options = ClaudeAgentOptions(
        system_prompt=orchestrator_prompt,
        allowed_tools=[
            "Read", "Write", "Edit", "Glob", "Grep",
            "Bash", "WebSearch", "WebFetch", "Task",
        ],
        permission_mode="acceptEdits",
        agents=agents,
        model="sonnet",
        cwd=str(config.default_cwd),
        mcp_servers={"skill-registry": skill_server},
    )

    print("Research Agent v0.2 — Domain Agents (type 'exit' to quit)\n")

    async with ClaudeSDKClient(options=options) as client:
        while True:
            try:
                user_input = input("You: ").strip()
            except (EOFError, KeyboardInterrupt):
                break

            if not user_input or user_input.lower() in ("exit", "quit"):
                break

            await client.query(user_input)
            async for msg in client.receive_response():
                handle_message(msg)
            print()


async def run_single(agent_name: str, prompt: str, config: AgentConfig) -> None:
    """Run a single domain agent directly."""
    agents = build_agents(config)

    if agent_name not in agents:
        available = ", ".join(sorted(agents))
        print(f"Unknown agent: {agent_name}\nAvailable: {available}")
        return

    agent_def = agents[agent_name]
    skill_server = create_skill_registry_server()

    options = ClaudeAgentOptions(
        system_prompt=agent_def.prompt,
        allowed_tools=agent_def.tools or [],
        permission_mode="acceptEdits",
        model=config.get_model(agent_name),
        cwd=str(config.default_cwd),
        mcp_servers={"skill-registry": skill_server},
    )

    async for msg in query(prompt=prompt, options=options):
        handle_message(msg)
