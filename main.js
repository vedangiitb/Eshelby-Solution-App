const { app, BrowserWindow } = require('electron');
const express = require('express');
const path = require('path');
const { spawn } = require('child_process');
const methodOverride = require('method-override');
const ejsMate = require('ejs-mate');
const bodyParser = require('body-parser');
const { appendFile } = require('fs');
const { error } = require('console');
// const { errors } = require('undici-types');

const appExpress = express();
const port = 3000;

//Body parser middleware
appExpress.use(bodyParser.urlencoded({extended:true}));

appExpress.engine('ejs',ejsMate);
appExpress.set('view engine', 'ejs');
appExpress.set('views', path.join(__dirname, 'views'));
appExpress.use(express.static(path.join(__dirname,'public')))
// appExpress.use(express.urlencoded({extended:true}))
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


appExpress.post('/isohomoinput',(req,res)=>{
  // const { a, b, c,eps11,eps22,eps33,eps12,eps13,eps23,ep,nu,mu } = req.body;
  const formData = req.body;
  // console.log(formData);

  const pythonProcess = spawn('python',['./Solution_codes/3D_isotropic_homogeneous.py', JSON.stringify(formData)]);
  let output = '';
  pythonProcess.stdout.on('data', (data) => {
    output += data.toString();
  });

  pythonProcess.on('close', (code) => {
    console.log(`Python script exited with code ${code}`);
    // Send the output back to the client
    console.log(output);
    res.render('result.ejs', { output });
    // res.send(output);
  });


  // try{
  //   // const pythonProcess = spawn('python',['./Solution_codes/3D_isotropic_homogeneous.py', a,b,c,eps11,eps22,eps33,eps12,eps13,eps23,ep,nu,mu])
  //   const pythonProcess = spawn('python',['./Solution_codes/3D_isotropic_homogeneous.py', JSON.stringify(formData)]);
  //   let output = " ";
  //   pythonProcess.stdout.on('data', (chunk) => {
  //     output += chunk.toString();
  //   // console.log(output)
  // });
  

  // pythonProcess.stderr.on('error', (error) => {
  //   console.error('Error from Python process:', error);
  //   res.status(500).send('Internal server error'); // Send an error response
  // });

  // // Once the Python process finishes, send the output back to the client
  // pythonProcess.on('close', (code) => {
  //   if (code === 0) {
  //       res.render('result.ejs', { output }); // Render the result page with processed data
  //   } else {
      
  //       console.error('Error code from Python process:', code);
  //       res.status(500).send('Internal server error');
  //   }
  // });
  // }

  // catch{
  //   console.error('Error handling form submission:', error);
  //   res.status(500).send('Internal server error');
  // }

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
