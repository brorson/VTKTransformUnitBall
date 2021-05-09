/*=========================================================================

  TransformUnitBall -- this calls vtkSphereSource to create & plot
  a unit ball (sphere) u.  It then creates a 3x3 random matrix and
  multiplies it into the unit ball, A*u, and plots the resulting
  ellipse.

  This code started from TestSphereSource.cxx, which is part of the
  VTK distribution.  I ripped out stuff I don't need and added
  my own suff.... lots of it.

  SDB -- March 2021.

=========================================================================*/

#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkAxesActor.h>
#include <vtkRenderWindow.h>
#include <vtkNamedColors.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkTextActor.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>
#include <vtkMatrix3x3.h>
#include <vtkCamera.h>

#include <iostream>
#include <ostream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector> 

//------------------------------------------------------------------
class KeyEventHandler {
public:
  bool RPressed;

  KeyEventHandler() {
    RPressed = false;
  }
  
  void KeypressCallbackFunction(vtkObject* caller,
                                long unsigned int vtkNotUsed(eventId),
                                void* vtkNotUsed(callData)) {
    vtkRenderWindowInteractor *iren = 
      static_cast<vtkRenderWindowInteractor*>(caller);

    //std::cout << "Caught event in KeypressCallbackFunction" << std::endl;
    std::string key;
    if (iren->GetKeySym()) {
      key = iren->GetKeySym();
      //cout << "Key pressed: " << key << endl;
    }
    
    if (key == "r"){
      RPressed = true;
      iren->ExitCallback ();  // Exit the loop since we caught an "r"
    } else {
      RPressed = false;
    }
  }

  bool RKeyPressed(void) {
    // I need to reset the state after every query, otherwise
    // the RPressed persists and I can't capture other events.
    bool temp = RPressed;
    RPressed = false;
    return temp;
  }
  
};

//=========================================================================
int main(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{

  cout << "Hit \'r\' to refresh (get a new ellipsoid)" << endl;
  cout << "Hit \'q\' to quit" << endl;
  
  // Seed RNG with time
  vtkMath::RandomSeed(time(NULL));

  //--------------------------------------------------------------------
  // Unit ball
  vtkSmartPointer<vtkSphereSource> unitBall = vtkSmartPointer<vtkSphereSource>::New();
  unitBall->SetThetaResolution(32);
  unitBall->SetPhiResolution(32);
  unitBall->SetStartTheta(0.0);
  unitBall->SetEndTheta(360.0);
  unitBall->SetStartPhi(0.0);
  unitBall->SetEndPhi(180.0);
  unitBall->LatLongTessellationOff(); // Generate triangles when off.
  unitBall->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);

  double center[3] = {0.0f, 0.0f, 0.0f};
  unitBall->SetCenter(center);

  double radius = 1.0f;
  unitBall->SetRadius(radius);

  unitBall->Update();

  //--------------------------------------------------------------------
  // Visualization stuff for unit ball
  
  // Define plot colors
  vtkNew<vtkNamedColors> colors;
  vtkColor3d actorColorRed = colors->GetColor3d("Red");
  vtkColor3d actorColorBlue = colors->GetColor3d("Blue");  
  vtkColor3d EdgeColour = colors->GetColor3d("slate_grey_dark");
  vtkColor3d BackgroundColour = colors->GetColor3d("Silver");
  vtkColor3d textColor = colors->GetColor3d("Green");  

  // Create mapper for unit ball
  vtkNew<vtkPolyDataMapper> unitBallMapper;
  unitBallMapper->SetInputConnection(unitBall->GetOutputPort());
  
  // Create actor for unit ball
  vtkNew<vtkActor> unitBallActor;
  unitBallActor->SetMapper(unitBallMapper);
  unitBallActor->GetProperty()->EdgeVisibilityOn();
  unitBallActor->GetProperty()->SetColor(actorColorBlue.GetData());
  unitBallActor->GetProperty()->SetEdgeColor(EdgeColour.GetData());
  unitBallActor->GetProperty()->SetOpacity(.3);  


  //-------------------------------------------------------------
  // General visualization stuff -- mostly boilerplate code.
  
  // Create renderer and render window
  vtkNew<vtkRenderer> renderer;
  renderer->SetBackground(BackgroundColour.GetData());

  // Store my initial camera position so I can get it back.
  double initPos[3] = {1.0, 0.0, 0.0};
  double initFP[3];
  double initViewUp[3] = {0.0, 0.0, 1.0}; // I want z axis up
  renderer->GetActiveCamera()->SetPosition(initPos);
  renderer->GetActiveCamera()->GetFocalPoint(initFP);
  renderer->GetActiveCamera()->SetViewUp(initViewUp);
  // cout << "Camera position = " << initPos[0] << " " << initPos[1] << " " << initPos[2] << " " << endl;
  
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Transforming a unit ball -- y = A*x");
  renderWindow->SetSize(800, 600); // This sets size of window on computer screen.
  
  
  // Create Interactor
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);
  vtkNew<vtkInteractorStyleTrackballCamera> style;
  renderWindowInteractor->SetInteractorStyle(style);
  
  // Key event -- attach to renderWindowInteractor
  KeyEventHandler KeyEvent;
  renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent,
				      &KeyEvent,
				      &KeyEventHandler::KeypressCallbackFunction);
  
  // Add axes to the render window interactor.
  // I do this new with every loop so the axes start in the
  // correct orientation.
  vtkNew<vtkAxesActor> axes;
  vtkNew<vtkOrientationMarkerWidget> widget;
  widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
  widget->SetOrientationMarker( axes );
  widget->SetInteractor( renderWindowInteractor );
  widget->SetViewport( 0.0, 0.0, 0.2, 0.2 );
  widget->SetEnabled( 1 );
    


  //-----------------------------------------------------------------
  // Now run a loop.  Inside the loop create a new random matrix
  // and transform the unit ball and make a plot of the resulting
  // ellipsoid.
  while(1) {

    cout << "----------------------------------------" << endl;
    //--------------------------------------------------------------------
    // Random matrix A.  Rather than pfaffing around with classes
    // provided by VTK I just used a standard C (not C++) 3x3 matrix
    // since it's easy.
    
    // Create random matrix.
    double val;
    double A[3][3];
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	val = vtkMath::Gaussian(0.0,1.0);
	A[i][j] = val;
      }
    }
    cout << "A = " << endl;
    for (int i=0; i<3; i++) {
      for (int j=0; j<2; j++) {
	cout << A[i][j] << "   ";
      }
      cout << A[i][2] << endl;
    }

    
    //-----------------------------------------------------------------
    // Compute SVD of A
    double U[3][3];
    double VT[3][3];
    double sigv[3];
    vtkMath::SingularValueDecomposition3x3 (A, U, sigv, VT);

    // Sort singular values in decreasing order.
    // Note that VTK uses a funky method to compute the SVD and
    // the singular values are not necesarily positive.  Therefore,
    // take the abs first.
    std::vector<double> sig(3);
    for (int i=0; i<3; i++) sig[i] = std::abs(sigv[i]);
    std::sort(sig.begin(), sig.end(), std::greater<double>());
    
    //cout << "sig = " << sig[0] << sig[1] << sig[2] << endl;
    std::stringstream out1;
    out1 << std::setprecision(3) << "sig1 = " << sig[0];
    std::stringstream out2;
    out2 << std::setprecision(3) << ", sig2 = " << sig[1];
    std::stringstream out3;
    out3 << std::setprecision(3) << ", sig3 = " << sig[2];
    //std::string txt1 = "sig1 = " + std::to_string(sig[0]);
    //std::string txt2 = "sig2 = " + std::to_string(sig[1]);
    //std::string txt3 = "sig3 = " + std::to_string(sig[2]);    
    //std::string txt = txt1 + ", " + txt2 + ", " +txt3;
    std::string txt = out1.str()+out2.str()+out3.str();
    cout << txt << endl;


    //--------------------------------------------------------------------
    // Create ellipsoid by multiplying unit ball by A.  I call
    // the result "ellipse".
    
    // Get points in unit ball, multiply by A and put result into new sphere object.
    vtkSmartPointer<vtkPoints> sphPts = vtkSmartPointer<vtkPoints>::New();
    sphPts = unitBall->GetOutput()->GetPoints();
    int N = sphPts->GetNumberOfPoints();
    //cout << "Number of points in unit ball N = " << N << endl;

    // Input and output vectors for transform.
    double x[3];  // Input vector to transformation
    double y[3];  // Output vector from transformation

    // Geometric pieces of the new ellipsoid.
    vtkSmartPointer<vtkPoints> ellipsePts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> ellipseVerts = vtkSmartPointer<vtkCellArray>::New();

    // The topology (mesh) of the ellipsoid is the same as the unit ball.  Use
    // the polygons in the unit ball as the polygons in the new ellipsoid.
    vtkSmartPointer<vtkCellArray> sphPolys = vtkSmartPointer<vtkCellArray>::New();  
    sphPolys = unitBall->GetOutput()->GetPolys();

    // Do the multiplication y = A*x
    double s;
    double ymax = 0;  // I use this to scale the actors when rendering.
    for (int i=0; i<N; i++) {
      sphPts->GetPoint(i, x);
      //cout << x[0] << ", " << x[1] << ", " << x[2] << endl;
      // Naive gemv
      for (int j=0; j<3; j++) {
	s = 0;
	for (int k=0; k<3; k++) {
	  s += A[j][k]*x[k];
	}
	y[j] = s;
	if (s > ymax) {
	  ymax = s;
	}
      }
      //cout << y[0] << ", " << y[1] << ", " << y[2] << endl;    
      //cout << "------------------" << endl;
      
      // Now put y into the new VTK object
      vtkIdType pid;    
      pid = ellipsePts->InsertNextPoint(y[0], y[1], y[2]);
      ellipseVerts->InsertNextCell( 1, &pid );
    }
    
    // Create a polydata object and add to it: points, verts, and polys.
    // If I don't add the polys then the plot only shows vertices, no
    // surface.
    vtkSmartPointer<vtkPolyData> ellipseData = 
      vtkSmartPointer<vtkPolyData>::New();
    ellipseData->SetPoints(ellipsePts);
    ellipseData->SetVerts (ellipseVerts);
    ellipseData->SetPolys(sphPolys);   // Reuse unit ball's polygon topology
    
    //-------------------------------------------------------------
    // Visualization stuff for ellipse
    
    // Create mapper for ellipse
    vtkSmartPointer<vtkPolyDataMapper> ellipseMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    ellipseMapper->SetInputData(ellipseData);
    ellipseMapper->ScalarVisibilityOn();
  
    // Create actor for ellipse
    vtkNew<vtkActor> ellipseActor;
    ellipseActor->SetMapper(ellipseMapper);
    ellipseActor->GetProperty()->EdgeVisibilityOn();
    ellipseActor->GetProperty()->SetColor(actorColorRed.GetData());
    ellipseActor->GetProperty()->SetEdgeColor(EdgeColour.GetData());
    ellipseActor->GetProperty()->SetOpacity(0.8);
    ellipseActor->GetProperty()->SetPointSize(1);

    //-------------------------------------------------------------
    // Text annotation
    vtkSmartPointer<vtkTextActor> textActor =
      vtkSmartPointer<vtkTextActor>::New();
    // txt is generated above -- the three sing. values.
    std::string legend = "Action of matrix on unit ball.  "+txt;
    textActor->SetInput(legend.c_str());
    textActor->GetProperty()->SetColor(textColor.GetData());
    textActor->SetPosition2(50, 50);
    textActor->GetTextProperty()->SetFontSize(24);
  
    // add the unit ball actor to the renderer
    double mag = 1.0/(2.2*(ymax+1.0));
    unitBallActor->SetScale(mag, mag, mag);
    renderer->AddActor(unitBallActor);

    // add the generated ellipse actor to the renderer
    ellipseActor->SetScale(mag, mag, mag);
    renderer->AddActor(ellipseActor);

    // Add the text to the renderer
    renderer->AddActor2D(textActor);

    // Start the event loop by telling the render window to render
    // and the interactor to start.  First reset view position.
    renderer->GetActiveCamera()->SetPosition(initPos);
    renderer->GetActiveCamera()->SetFocalPoint(initFP);
    renderer->GetActiveCamera()->SetViewUp(initViewUp);
    renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();
    
    // If we get here, it's because an event happened which caused
    // the interactor to stop.  Check if it's an 'r' key for refresh
    // or if some other even happened, in case we should quit.
    if(KeyEvent.RKeyPressed()) {
      //cout << "R key noticed in main" << endl;
      renderer->RemoveActor(ellipseActor);
      renderer->RemoveActor2D(textActor);
    } else {
      //cout << "Back in main, no R key ... must be time to quit." << endl;
      return EXIT_SUCCESS;
    }

  } // while(1)
    
}
