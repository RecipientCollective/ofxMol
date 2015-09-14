#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    ofBackground(0);
    
    // setup system
    system.setup(ofFile("ice.pdb").path());
    mesh = system.getModel(0).atomsMesh(0.5f);
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
    ofTranslate(ofGetWidth()/2.0, ofGetHeight()/2.0);
    ofScale(40.0, 40.0);
    ofRotate(ofGetElapsedTimef()*1.5,0,1,1);
    
    // draw atoms
    ofSetColor(255);
    mesh.draw();
    ofPopMatrix();

    ofDisableDepthTest();
    
    // log info
    ofSetColor(255);
    ofDrawBitmapString(system.getModel(0).log(),10,20);
}