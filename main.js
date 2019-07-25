const {app, BrowserWindow} = require('electron');

function createWindow () {

    let python = require('child_process').spawn('python', ['./__init__.py']);
    let window = new BrowserWindow({width: 1366, height: 768});
    // window.loadURL('index.html');
    window.loadURL('http://127.0.0.1:5000/');
}
app.on('ready', createWindow);
app.on('window-all-closed', () => {
// On macOS it is common for applications and their menu bar
// to stay active until the user quits explicitly with Cmd + Q
if (process.platform !== 'darwin') {
  app.quit()
}
});
