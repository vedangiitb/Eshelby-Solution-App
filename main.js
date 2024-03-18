const { app, BrowserWindow, ipcMain} = require('electron');
const express = require('express');
const path = require('path');
const { spawn } = require('child_process');
const methodOverride = require('method-override');
const ejsMate = require('ejs-mate');
const bodyParser = require('body-parser');

const appExpress = express();
const port = 3000;

//Body parser middleware
appExpress.use(bodyParser.urlencoded({extended:true}));

appExpress.engine('ejs',ejsMate);
appExpress.set('view engine', 'ejs');
appExpress.set('views', path.join(__dirname, 'views'));
appExpress.use(express.static(path.join(__dirname,'public')))
appExpress.use(methodOverride('_method'));

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

appExpress.get('/learnproblem',(req,res)=>{
  res.render('Tuts/problemTut')
})

appExpress.get('/learnsoftware',(req,res)=>{
  res.render('Tuts/softwareTut')
})

appExpress.get('/settings',(req,res)=>{
  res.render('settings')
})

appExpress.post('/isohomoinput',(req,res)=>{
  const formData = req.body;

  const pythonProcess = spawn('python',['./Solution_codes/test.py', JSON.stringify(formData)]);
  let output = '';
  pythonProcess.stdout.on('data', (data) => {
    output += data.toString();
  });

  pythonProcess.on('close', (code) => {
    console.log(`Python script exited with code ${code}`);
    console.log(output);
    res.render('result.ejs', { output });
  });
})

const server = appExpress.listen(port, () => {
  console.log(`Express server running at http://localhost:${port}`);
});

function createWindow() {
  const win = new BrowserWindow({
    autoHideMenuBar: true,
    width: 850,
    height: 600,
    webPreferences: {
      nodeIntegration: true
    }
  });

  win.loadURL(`http://localhost:${port}`);

  ipcMain.on('goBack', () => {
    if (mainWindow && mainWindow.webContents.canGoBack()) {
      mainWindow.webContents.goBack();
    }
  });
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
