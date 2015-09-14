#pragma once

#include "ofMain.h"
#include "ofxMol.h"

class ofApp : public ofBaseApp
{
    
public:
    void setup();
    void update();
    void draw();
    
protected:
    OfxMol::System system;
    ofMesh atoms;
};
