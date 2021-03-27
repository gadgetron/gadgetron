# Dev Container for VS Code

The `.devcontainer` configuration enabled [development in a container](https://code.visualstudio.com/docs/remote/containers).

After starting in the container, post create commands will attempt to populate your `.vscode` folder. You can also copy `.devcontainter/.vscode` to your workspace:

```bash
.devcontainer/bootstrap-vscode.sh 
```

This will give you a starting set of VS code tasks. You should then be able to simply hit F5 to build and start in debug mode.