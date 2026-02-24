# Project Guidelines

## Agents

### /review-bloat — Code Bloat Review Agent

Run with: `/review-bloat` or `/review-bloat path/to/file.py`

Scan the codebase (or a specific file) for bloated, overly complex, or hard-to-read code. The goal is simple, readable code.

**What to flag:**

- Dead code: unused imports, unreachable branches, commented-out code, variables assigned but never read
- Over-abstraction: wrapper functions that add no logic, single-use helper classes, unnecessary inheritance hierarchies
- Redundant logic: repeated null/type checks already guaranteed by callers, defensive code for impossible states, duplicate validation at multiple layers
- Verbose patterns: manually building something the stdlib already provides, reimplementing builtins (e.g., hand-rolled CSV parsing when `csv` module works), long if/elif chains replaceable by a dict lookup
- Complex control flow: deeply nested conditionals (3+ levels), functions longer than ~50 lines that do multiple things, boolean parameters that switch behavior (should be separate functions)
- Unnecessary flexibility: config options nobody changes, parameterized code with only one call site, premature generalization

**How to report:**

For each finding, output:
1. File and line range
2. What the issue is (one sentence)
3. A concrete fix (refactored code snippet or deletion)

Prioritize by impact: large deletions and simplifications first, style nits last. Skip test files unless explicitly asked. Do not flag docstrings, type annotations, or logging unless they are clearly excessive.
