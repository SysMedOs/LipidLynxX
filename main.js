const {app, BrowserWindow} = require('electron');

python_worker = require('child_process').spawn('python', ['./lynx/app.py']);
console.log(`Start Python worker...`);

function createWindow () {
    let window = new BrowserWindow({width: 1366, height: 800});
    window.setMenuBarVisibility(false);
    // load the local epiLION website powered by flask
    window.loadURL('http://127.0.0.1:5000/');
    console.log(`local service ready...`);
}
function startLION() {
  console.log(`Start epiLION...`);
  // Wait few seconds to let python flask web site fully functional
    // before starting main window
  setTimeout(createWindow, 3333);
}

app.on('ready', startLION);
app.on('window-all-closed', () => {
    app.quit();
  // make sure that python & flask closed properly
  python_worker.kill();
});
