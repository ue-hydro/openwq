AI-Powered Assistance
=====================

OpenWQ provides built-in support for AI-powered assistance using Claude and other large language models. Two specialized prompt files help AI assistants understand the OpenWQ codebase and provide accurate, context-aware help.

Getting Started
~~~~~~~~~~~~~~~

**For Claude Code CLI Users (Users & Developers)**

Claude Code is a command-line AI assistant that can read, write, and understand code directly in your terminal. The ``CLAUDE.md`` file in the OpenWQ repository root is automatically loaded, giving Claude full context about the project.

**Installation (macOS/Linux):**

.. code-block:: bash

    # Step 1: Download and install Claude Code
    curl -fsSL https://claude.ai/install.sh | sh

    # Step 2: Add Claude to your PATH (REQUIRED - copy and run this entire command)
    # For zsh (default on macOS):
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc && source ~/.zshrc

    # For bash (common on Linux):
    # echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc

    # Step 3: Verify installation works
    claude --version

    # Step 4: Log in with your Claude account (first time only)
    claude

**Installation (Windows):**

.. code-block:: powershell

    # Step 1: Install using winget (recommended)
    winget install Anthropic.ClaudeCode

    # Alternative: Download installer from https://github.com/anthropics/claude-code/releases

    # Step 2: Close and reopen your terminal (required for PATH to update)

    # Step 3: Verify installation works
    claude --version

    # Step 4: Log in with your Claude account (first time only)
    claude

**Troubleshooting:**

If ``claude`` command is not found after installation:

* **macOS/Linux:** Run the PATH command from Step 2 above, then try again
* **Windows:** Make sure you completely closed and reopened your terminal
* **All platforms:** Try opening a brand new terminal window

**Usage:**

.. code-block:: bash

    # Navigate to OpenWQ directory
    cd /path/to/openwq/

    # Start Claude Code in interactive mode
    claude

    # Claude automatically reads CLAUDE.md and understands the project
    # Ask questions directly in the terminal:
    > How do I set up a nitrification reaction?
    > Help me configure the calibration framework
    > What files handle sediment transport?

    # Or run a single command without interactive mode:
    claude "Explain how the DDS optimizer works"

**What Claude Code Can Do:**

- Read and understand any file in the repository
- Write and edit code files directly
- Run shell commands (with your permission)
- Search across the codebase
- Explain complex code sections
- Debug errors and suggest fixes

**Example Session:**

.. code-block:: bash

    $ cd openwq/
    $ claude

    Claude Code v1.x.x
    Working directory: /path/to/openwq
    Loaded: CLAUDE.md (project context)

    You: How do I add a new BGC reaction for phosphorus uptake?

    Claude: I'll help you add a phosphorus uptake reaction. Based on the
    NATIVE_BGC_FLEX templates in supporting_scripts/Model_Config/...
    [Claude provides detailed instructions and can edit files directly]

**For Web/API Users (All Users)**

Copy the prompt from ``docs/OPENWQ_ASSISTANT_PROMPT.md`` into your AI conversation:

1. Open ``docs/OPENWQ_ASSISTANT_PROMPT.md``
2. Copy the entire content
3. Paste as your first message in Claude (claude.ai) or other AI assistants
4. Then ask your OpenWQ questions


Claude Project Setup (Recommended for Teams)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Claude Pro/Team subscribers, create a dedicated project:

1. Go to claude.ai → Projects → Create Project
2. Name it "OpenWQ Assistant"
3. In Project Instructions, paste the content from ``docs/OPENWQ_ASSISTANT_PROMPT.md``
4. Upload key documentation files:

   - ``CLAUDE.md``
   - ``calibration_config_template.py``
   - Key JSON configuration examples
   - Relevant RST documentation files

5. All conversations in this project will have full OpenWQ context


What the AI Assistant Can Help With
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Configuration & Setup**

- Writing BGC reaction configurations (NATIVE_BGC_FLEX, PHREEQC)
- Setting up sorption isotherms (Freundlich, Langmuir)
- Configuring sediment transport modules
- Creating source/sink load files

**Calibration**

- Setting up the calibration framework
- Choosing objective functions and temporal resolution
- Configuring sensitivity analysis (Morris, Sobol)
- Interpreting calibration results

**Troubleshooting**

- Debugging configuration errors
- Understanding error messages
- Fixing compilation issues
- Resolving HPC/container problems

**Code Understanding**

- Explaining C++ source code structure
- Understanding Python supporting scripts
- Navigating the repository


File Locations
~~~~~~~~~~~~~~

+----------------------------------------+------------------------------------------------+
| File                                   | Purpose                                        |
+========================================+================================================+
| ``CLAUDE.md``                          | Auto-loaded by Claude Code CLI                 |
+----------------------------------------+------------------------------------------------+
| ``docs/OPENWQ_ASSISTANT_PROMPT.md``    | Copy-paste prompt for web/API users            |
+----------------------------------------+------------------------------------------------+


API Integration
~~~~~~~~~~~~~~~

For programmatic use with the Anthropic API:

.. code-block:: python

    import anthropic

    # Load the OpenWQ system prompt
    with open("docs/OPENWQ_ASSISTANT_PROMPT.md") as f:
        system_prompt = f.read()

    client = anthropic.Anthropic()
    message = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=4096,
        system=system_prompt,
        messages=[
            {"role": "user", "content": "How do I configure denitrification?"}
        ]
    )
    print(message.content[0].text)


Tips for Best Results
~~~~~~~~~~~~~~~~~~~~~

1. **Be specific** - Include file names, species names, and error messages
2. **Share context** - Paste relevant JSON snippets or code sections
3. **Iterate** - Ask follow-up questions to refine solutions
4. **Verify** - Always validate AI suggestions against documentation


Updating the Prompts
~~~~~~~~~~~~~~~~~~~~

The AI assistant prompts should be updated when:

- New modules are added to OpenWQ
- Configuration file formats change
- New supporting scripts are created
- Common troubleshooting patterns emerge

To update, edit ``CLAUDE.md`` and ``docs/OPENWQ_ASSISTANT_PROMPT.md`` and commit the changes.
