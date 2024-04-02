const { app, BrowserWindow, ipcMain} = require('electron');
const express = require('express');
const path = require('path');
const { spawn } = require('child_process');
const methodOverride = require('method-override');
const ejsMate = require('ejs-mate');
const bodyParser = require('body-parser');
const appExpress = express();
const port = 3000;


const validateInput = (req, res, next) => {
  const { a, b, c } = req.body;
  const aValue = parseFloat(a);
  const bValue = parseFloat(b);
  const cValue = parseFloat(c);

  if (aValue > bValue && bValue > cValue) {
    next();
  } else {
    res.send("Invalid Input");
  }
};
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

appExpress.post('/isohomoinput',validateInput,(req,res)=>{
  const formData = req.body;

  const pythonProcess = spawn('python',['./Solution_codes/3D_isotropic_homogeneous.py', JSON.stringify(formData)]);
  let output = '';
  let error = '';
  pythonProcess.stdout.on('data', (data) => {
    output += data.toString();
  });

      
  pythonProcess.stderr.on('data', (data) => {
    error += data.toString();
  });

  pythonProcess.on('close', (code) => {
    console.log(`Python script exited with code ${code}`);
    console.log(error);
    console.log(output);
    res.render('result.ejs', { output });
  });

})

appExpress.get('/dummy',(req,res)=>{
  const pythonProcess = spawn('python',['./Solution_codes/dummy.py','plot']);

  // Handle Python process events
  pythonProcess.stdout.on('data', (data) => {
    console.log(`stdout: ${data}`);
  });

  pythonProcess.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
  });

  pythonProcess.on('close', (code) => {
    console.log(`Python process exited with code ${code}`);
    res.send('Python function executed successfully');
  });
});

appExpress.post('/saveData', (req, res) => {
  const data = req.body.data; // Assuming the data is sent in the request body
  const fileName = 'savedData.txt';

  // Write data to a text file
  fs.writeFile(fileName, data, (err) => {
    if (err) {
      console.error(err);
      res.status(500).send('Error saving data');
    } else {
      console.log('Data saved successfully');
      res.send('Data saved successfully');
    }
  });
});


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
