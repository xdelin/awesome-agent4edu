# Support Boundaries

This document sets expectations for issue triage, support scope, and rollout safety.

## In Scope

- Reproducible bugs in published release behavior.
- Install and setup blockers for supported Node/runtime paths.
- Regressions in existing tools (`slack_get_thread`, `slack_search_messages`, etc.).
- Documentation corrections with clear technical evidence.

## Out of Scope

- Custom Slack policy/legal interpretation for a specific organization.
- Workspace-specific data migrations or bespoke integrations.
- Real-time operational support for third-party hosting providers.
- Urgent production incident ownership for self-hosted external deployments.

## Intake Requirements

Include the following in every issue:

1. Version (`npm view @jtalk22/slack-mcp version` output).
2. Runtime mode (`stdio`, `web`, `http`, or Worker/Smithery).
3. OS and Node version.
4. Minimal reproduction steps and exact error text.
5. Whether this blocks individual use or team rollout.

## Response Windows (Best Effort)

- Installation/blocker bugs: initial response target within 2 business days.
- Non-blocking bugs: initial response target within 5 business days.
- Feature requests: triaged in backlog batches.

## Solo Maintainer Capacity Guardrail

- Weekly support budget target: <= 2 hours/week.
- If inbound support exceeds this cap, new feature work may be paused.
- High-context requests may be deferred until reproducible artifacts are provided.

## Security and Data Handling

- Do not post raw tokens/cookies in issues.
- Use redacted logs.
- If credentials are exposed, rotate immediately and update issue with redaction.

## Deployment Escalation Rule

For team/hosted usage, open a deployment intake issue before broad rollout.
