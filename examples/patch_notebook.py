import json
from pathlib import Path

notebook_path = Path("specviz_full_pipeline.ipynb")
with open(notebook_path, "r") as f:
    nb = json.load(f)

# The code to find and replace
target_start = '# Try to find specviz directory\n'
replacement_code = [
    '# Try to find specviz directory\n',
    'if (cwd / "app.py").exists() and (cwd / "utils").exists():\n',
    '    specviz_dir = cwd\n',
    'elif (cwd.parent / "app.py").exists() and (cwd.parent / "utils").exists():\n',
    '    specviz_dir = cwd.parent\n',
    'else:\n',
    '    # Search parents for project structure (app.py + utils folder)\n',
    '    for parent in [cwd] + list(cwd.parents):\n',
    '        if (parent / "app.py").exists() and (parent / "utils").exists():\n',
    '            specviz_dir = parent\n',
    '            break\n'
]

# Find the cell
modified = False
for cell in nb["cells"]:
    if cell["cell_type"] == "code":
        source = cell["source"]
        # Look for the block
        try:
            start_idx = -1
            end_idx = -1
            for i, line in enumerate(source):
                if line.strip() == '# Try to find specviz directory':
                    start_idx = i
                if start_idx != -1 and line.strip() == 'break':
                    # This is a bit brittle, relying on the break at the end of the else block
                    # Let's check the context from the file view
                    # The original ending was:
                    #             specviz_dir = parent / "specviz"
                    #             break
                    pass
            
            # Since the replacement is a clean block replacement, let's just find the start and identifying end
            # The original code block:
            # 112: if (cwd / "specviz").exists():
            # ...
            # 122:             break
            
            # Let's find the index of the start line
            for i, line in enumerate(source):
                if line.strip() == '# Try to find specviz directory':
                    start_idx = i
                    break
            
            if start_idx != -1:
                # Found the start.
                # The block ends before "if specviz_dir is None:"
                end_idx = -1
                for i in range(start_idx, len(source)):
                    if "if specviz_dir is None:" in source[i]:
                        end_idx = i
                        break
                
                if end_idx != -1:
                    # Remove the old lines (excluding the next block start)
                    # We want to replace from start_idx to end_idx (exclusive)
                    # Use a slice assignment
                    cell["source"][start_idx:end_idx] = replacement_code + ['\n']
                    modified = True
                    print("✅ Found and patched the cell.")
                    break
        except Exception as e:
            print(f"Error processing cell: {e}")

if modified:
    with open(notebook_path, "w") as f:
        json.dump(nb, f, indent=1)
    print("✅ Notebook saved.")
else:
    print("❌ Could not find the target code block to patch.")
