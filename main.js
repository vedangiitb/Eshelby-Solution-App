const { app, BrowserWindow, ipcMain,shell} = require('electron');
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

let inputData = {};

function convertToMM(value, units) {
  switch (units) {
      case 'cm':
          return value * 10; 
      case 'm':
          return value * 1000; 
      default:
          return value; // No conversion needed for default unit (mm)
  }
}

function convertToGPa(value, units) {
  switch (units) {
      case 'MPa':
          return value * 0.001; 
      default:
          return value; // No conversion needed for default unit(GPa)
  }
}

appExpress.post('/isohomoinput',validateInput,(req,res)=>{

  // Extract non-dynamic input
  const a = req.body.a;
  const b = req.body.b;
  const c = req.body.c;
  const length_units = req.body.length_units;
  const eps11 = req.body.eps11;
  const eps22 = req.body.eps22;
  const eps33 = req.body.eps33;
  const eps13 = req.body.eps13;
  const eps23 = req.body.eps23;
  const eps12 = req.body.eps12;
  const ep = req.body.ep;
  const E_units = req.body.E_units;
  const nu = req.body.nu;



  // Extract dynamic inputs

  const inputCount = parseInt(req.body.inputCount);

  const dynamicInputs = [];
    for (let i = 1; i <= inputCount; i++) {
        const values = [];
        
            const x = req.body[`value${i}_x`];
            const y = req.body[`value${i}_y`];
            const z = req.body[`value${i}_z`];
            
            values.push(convertToMM(parseFloat(x), length_units));
            values.push(convertToMM(parseFloat(y), length_units));
            values.push(convertToMM(parseFloat(z), length_units));
        
        dynamicInputs.push(values);
    }
  

  // Convert values to millimeters
  const a_mm = convertToMM(parseFloat(a),length_units );
  const b_mm = convertToMM(parseFloat(b), length_units);
  const c_mm = convertToMM(parseFloat(c), length_units);
  const E_GPa = convertToGPa(parseFloat(ep), E_units);

  // Update inputData
  inputData = {
    a: a_mm,
    b: b_mm,
    c: c_mm,
    eps11:eps11,
    eps22:eps22,
    eps33:eps33,
    eps13:eps13,
    eps23:eps23,
    eps12:eps12,
    ep:E_GPa,
    nu:nu,
    targets:dynamicInputs

};

  console.log(inputData);


  const pythonProcess = spawn('python',['./Solution_codes/3D_isotropic_homogeneous.py', JSON.stringify(inputData)]);
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

appExpress.post('/isoinhomoinput',validateInput,(req,res)=>{

  // Extract non-dynamic input
  const a = req.body.a;
  const b = req.body.b;
  const c = req.body.c;
  const length_units = req.body.length_units;
  const eps11 = req.body.eps11;
  const eps22 = req.body.eps22;
  const eps33 = req.body.eps33;
  const eps13 = req.body.eps13;
  const eps23 = req.body.eps23;
  const eps12 = req.body.eps12;
  const ep = req.body.ep;
  const eps = req.body.eps;
  const E_units = req.body.E_units;
  const nu = req.body.nu;
  const nus = req.bodu.nus;

  // Extract dynamic inputs
  const inputCount = parseInt(req.body.inputCount);

  const dynamicInputs = [];
    for (let i = 1; i <= inputCount; i++) {
        const values = [];
        
            const x = req.body[`value${i}_x`];
            const y = req.body[`value${i}_y`];
            const z = req.body[`value${i}_z`];
            
            values.push(convertToMM(parseFloat(x), length_units));
            values.push(convertToMM(parseFloat(y), length_units));
            values.push(convertToMM(parseFloat(z), length_units));
        
        dynamicInputs.push(values);
    }
  

  // Convert values to millimeters
  const a_mm = convertToMM(parseFloat(a),length_units );
  const b_mm = convertToMM(parseFloat(b), length_units);
  const c_mm = convertToMM(parseFloat(c), length_units);
  const E_GPa = convertToGPa(parseFloat(ep), E_units);
  const ES_GPa = convertToGPa(parseFloat(eps),E_units);

  // Update inputData
  inputData = {
    a: a_mm,
    b: b_mm,
    c: c_mm,
    eps11:eps11,
    eps22:eps22,
    eps33:eps33,
    eps13:eps13,
    eps23:eps23,
    eps12:eps12,
    ep:E_GPa,
    eps: ES_GPa,
    nu:nu,
    nus:nus,
    targets:dynamicInputs
};

  console.log(inputData);


  const pythonProcess = spawn('python',['./Solution_codes/3D_isotropic_inhomogeneous.py', JSON.stringify(inputData)]);
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

  win.maximize();

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
