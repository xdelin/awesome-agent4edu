# test-generator

Universal test generation skill for Claude Code that learns from your project.

## Features

- 🔍 **Auto-detects** your language and test framework
- 🧠 **Self-learning** — improves from your modifications
- 🌍 **9 languages** — TypeScript, JavaScript, Python, PHP, Go, Java, C#, Ruby, Rust
- 📁 **Zero config** — just ask for tests

## Installation
```bash
# In Claude Code
/skill install https://github.com/mkdirrobert/test-generator
```

Or download `test-generator.skill` from [Releases](https://github.com/mkdirrobert/test-generator/releases) and:
```bash
/skill install /path/to/test-generator.skill
```

## Usage

Just ask naturally:

- "Write tests for UserService"
- "Generate unit tests for src/utils/validation.ts"
- "Create API tests for my login endpoint"

## How It Works

1. Detects your project's language from manifest files
2. Assesses quality of existing tests
3. Generates tests matching YOUR conventions
4. Learns from your corrections

## License

MIT
