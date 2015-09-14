#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    bShowHelp = true;
    bCoarseAtomsMesh = false;
    bAtomsPointCloud = true;
    bAtoms = false;
    bBackbone = false;
    
    // LOAD PDB FILES IN DATA DIR
    dataDir.open(DATADIR);
    ofLogNotice() << "FILES:" << dataDir.listDir();
    currentFile = 0;
    
    if (!dataDir.exists())
    {
        ofLogFatalError() << "The directory: " << dataDir.path() << " is missing";
        ofExit();
    }
    
    dataDir.allowExt("pdb");
    pdbFiles = dataDir.getFiles();
    
    for (std::vector<ofFile>::iterator fit=pdbFiles.begin(); fit!=pdbFiles.end(); ++fit)
    {
        ofLogNotice() << "PDB FILE: " << fit->getFileName();
    }

    
    ofSetFrameRate(60);
    
    ofSetVerticalSync(true);
    cam.setDistance(100);
    ofSetCircleResolution(64);
    ofBackground(0);
    
    pointLightLeft.setPosition(-100, 20, 0);
    pointLightRight.setPosition(100, -20, 0);
    
    pointLightLeft.setDiffuseColor( ofFloatColor(.35, .35, .85) );
    pointLightLeft.setSpecularColor( ofFloatColor(1.f, 1.f, 1.f));
    
    pointLightRight.setDiffuseColor( ofFloatColor(.35, .35, .85) );
    pointLightRight.setSpecularColor( ofFloatColor(1.f, 1.f, 1.f));

    loadMolecule(pdbFiles[currentFile]);
}

//--------------------------------------------------------------
void ofApp::update()
{
    ofSetWindowTitle(ofToString(pdbFiles[currentFile].getFileName()));
}

void ofApp::loadMolecule(ofFile file)
{
    system.setup(file.path(), OfxMol::ADVANCED);
    cAtmMesh = system.getModel(0).coarseAtomsMesh(8);
    atoms = system.getModel(0).atomsMesh(1.0f);
    atomsCloud = system.getModel(0).atomsPointCloud();
    backbone = system.getModel(0).backbonePoly();
}

//--------------------------------------------------------------
void ofApp::draw()
{
    ofEnableDepthTest();
    cam.begin();
    
    ofRotateX(ofRadToDeg(.5));
    ofRotateY(ofRadToDeg(-.5));
    
    ofEnableLighting();

    pointLightLeft.enable();
    pointLightRight.enable();
    
    material.begin();
    ofSetColor(255);
    
    ofPushMatrix();
    ofTranslate(0.0f, 0.0f,0.0f);
    
    if (bCoarseAtomsMesh)
    {
        cAtmMesh.draw();
    }
    
    if (bAtomsPointCloud)
    {
        glPointSize(24.0);
        atomsCloud.drawVertices();
    }
    
    if (bAtoms)
    {
        atoms.draw();
    }
    
    if (bBackbone)
    {
        glLineWidth(4.0f);
        backbone.draw();
    }
    
    ofPopMatrix();
    
    material.end();
    
    cam.end();
    
    ofDisableLighting();
    ofDisableDepthTest();
    
    drawInteractionArea();
    ofSetColor(255);
    string msg = string("Using mouse inputs to navigate (press 'c' to toggle): ") + (cam.getMouseInputEnabled() ? "YES" : "NO");
    msg += string("\nShowing help (press 'h' to toggle): ")+ (bShowHelp ? "YES" : "NO");
    
    if (bShowHelp)
    {
        msg += "\n\nATOMS VISUALIZATION:\n";
        msg += "Coarse Atoms Mesh [1] : " + ofToString(bCoarseAtomsMesh ? "YES" : "NO") + "\n";
        msg += "atoms point cloud [2] : " + ofToString(bAtomsPointCloud ? "YES" : "NO") + "\n";
        msg += "atoms spheres [3] : " + ofToString(bAtoms ? "YES" : "NO") + "\n";
        msg += "backbone [4] : " + ofToString(bBackbone ? "YES" : "NO") + "\n";
        msg += "\n\nLEFT MOUSE BUTTON DRAG:\nStart dragging INSIDE the yellow circle -> camera XY rotation .\nStart dragging OUTSIDE the yellow circle -> camera Z rotation (roll).\n\n";
        msg += "LEFT MOUSE BUTTON DRAG + TRANSLATION KEY (" + ofToString(cam.getTranslationKey()) + ") PRESSED\n";
        msg += "OR MIDDLE MOUSE BUTTON (if available):\n";
        msg += "move over XY axes (truck and boom).\n\n";
        msg += "RIGHT MOUSE BUTTON:\n";
        msg += "move over Z axis (dolly)";
    }
    
    msg += "\n\nfps: " + ofToString(ofGetFrameRate(), 2);
    ofDrawBitmapStringHighlight(msg, 10, 20);

}

//--------------------------------------------------------------
void ofApp::drawInteractionArea()
{
    ofRectangle vp = ofGetCurrentViewport();
    float r = MIN(vp.width, vp.height) * 0.5f;
    float x = vp.width * 0.5f;
    float y = vp.height * 0.5f;
    
    ofPushStyle();
    ofSetLineWidth(3);
    ofSetColor(255, 255, 0);
    ofNoFill();
    glDepthMask(false);
    ofCircle(x, y, r);
    glDepthMask(true);
    ofPopStyle();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{
    switch(key)
    {
        case 'C':
        case 'c':
            if(cam.getMouseInputEnabled()) cam.disableMouseInput();
            else cam.enableMouseInput();
            break;
            
        case 'F':
        case 'f':
            ofToggleFullscreen();
            break;
        case 'H':
        case 'h':
            bShowHelp ^=true;
            break;
            
        case '1':
            bCoarseAtomsMesh ^=true;
            break;
        case '2':
            bAtomsPointCloud ^=true;
            break;
        case '3':
            bAtoms ^=true;
            break;
        case '4':
            bBackbone ^=true;
            break;
            
        case OF_KEY_LEFT:
            currentFile--;
            if (currentFile < 0)
            {
                currentFile = 0;
            }
            loadMolecule(pdbFiles[currentFile]);
            break;
        case OF_KEY_RIGHT:
            currentFile++;
            if (currentFile > pdbFiles.size()-1)
            {
                currentFile = pdbFiles.size()-1;
            }
            loadMolecule(pdbFiles[currentFile]);
            break;
    }
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
