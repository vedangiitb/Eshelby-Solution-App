const { app, BrowserWindow } = require('electron');
const express = require('express');
const path = require('path');

const appExpress = express();
const port = 3000;

appExpress.set('view engine', 'ejs');
appExpress.set('views', path.join(__dirname, 'views'));
appExpress.use(express.static(path.join(__dirname,'public')))

appExpress.get('/', (req, res) => {
  res.render('index');
});

appExpress.get('/type',(req,res)=>{
  res.render('problemType')
})

appExpress.get('/isohomogen',(req,res)=>{
  res.render('InputPages/isohomogen')
})

appExpress.get('/isohetergen',(req,res)=>{
  res.render('InputPages/isohetergen')
})

appExpress.get('/anisohomogen',(req,res)=>{
  res.render('InputPages/anisohomogen')
})

appExpress.get('/anisoheterogen',(req,res)=>{
  res.render('InputPages/anisoheterogen')
})

const server = appExpress.listen(port, () => {
  console.log(`Express server running at http://localhost:${port}`);
});

function createWindow() {
  const win = new BrowserWindow({
    width: 800,
    height: 600,
    webPreferences: {
      nodeIntegration: true
    }
  });

  win.loadURL(`http://localhost:${port}`);
}

app.whenReady().then(createWindow);

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createWindow();
  }
});
