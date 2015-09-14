#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    system.setup(ofFile("2WY4.pdb").path());
    ofBackground(0);
    line = system.getModel(0).backbonePoly();
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
    line.draw();
    ofPopMatrix();
    ofDisableDepthTest();
    
    ofSetColor(255, 255, 255);
    ofDrawBitmapString(system.getModel(0).log(),10,20);
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

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
