#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    system.setup(ofFile("2WY4.pdb").path(), OfxMol::ADVANCED);
    pointLight.setPosition(20.0f, 20.0f, 0.0f);
    pointLightRetro.setPosition(ofGetWidth() - 20.0f, ofGetHeight() - 20.0f, 0.0f);
    mesh = system.getModel(0).coarseAtomsMesh(ofColor(255), 12);
    ofLogNotice() << "MESH SIZE: " << mesh.getNumVertices() << " NORMALS: " << mesh.getNumNormals();
    ofLogNotice() << "COARSE ATOMS SIZE: " << system.getModel(0).number_of_coarse_atoms();
    
    ofSetSmoothLighting(true);
    pointLight.setDiffuseColor( ofFloatColor(.85, .85, .55) );
    pointLight.setSpecularColor( ofFloatColor(1.f, 1.f, 1.f));
    pointLightRetro.setDiffuseColor( ofFloatColor(.85, .85, .55) );
    pointLightRetro.setSpecularColor( ofFloatColor(1.f, 1.f, 1.f));
    
    // shininess is a value between 0 - 128, 128 being the most shiny //
    material.setShininess( 80 );
    // the light highlight of the material //
    material.setSpecularColor(ofColor(155, 155, 155, 255));
}

//--------------------------------------------------------------
void ofApp::update()
{
    ofSetWindowTitle(ofToString(ofGetFrameRate(), 0));
}

//--------------------------------------------------------------
void ofApp::draw()
{
    ofBackgroundGradient(ofColor(150),ofColor(100), OF_GRADIENT_BAR);
    
    ofEnableDepthTest();
    ofEnableLighting();
    pointLight.enable();
    pointLightRetro.enable();
    
    material.begin();
    
    ofPushMatrix();
    ofTranslate(ofGetWidth()/2.0f, ofGetHeight()/2.0f);
    ofScale(10,10);
    ofRotate(ofGetElapsedTimef()*10, 0, 1, 0);
    ofSetColor(100);
    mesh.draw();
    ofPopMatrix();
    
    material.end();
    
    ofDisableLighting();
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
