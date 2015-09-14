#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    system.setup(ofFile("2WY4.pdb").path());
    atoms = system.getModel(0).atomsPointCloud();
    ofBackground(0);
    glPointSize(4.0f);
    
    // log atom info for first 10 atoms
    for (int i = 0; i < 10; i++)
    {
        OfxMol::Atom atom = system.getModel(0).getAtom(i);
        ofLogNotice() << atom.log();
    }
}

//--------------------------------------------------------------
void ofApp::update()
{
    ofSetWindowTitle(ofToString(ofGetFrameRate(), 0));
}

//--------------------------------------------------------------
void ofApp::draw()
{
    ofEnableDepthTest();
    ofPushMatrix();
    ofTranslate(ofGetWidth()/2.0f, ofGetHeight()/2.0f);
    ofScale(10,10);
    ofRotate(ofGetElapsedTimef()*10, 0, 1, 0);
    atoms.drawVertices();
    ofPopMatrix();
    ofDisableDepthTest();
    
    ofSetColor(255, 255, 255);
    ofDrawBitmapString(system.getModel(0).log(),10,20);
}