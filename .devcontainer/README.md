# Dev Container for VS Code

The `.devcontainer` configuration enabled [development in a container](https://code.visualstudio.com/docs/remote/containers).

After starting in the container, you can copy `.devcontainter/.vscode` to your workspace:

```bash
cp -r .devcontainer/.vscode ./
```

This will give you a starting set of VS code tasks. You should then be able to simply hit F5 to build and start in debug mode.