#pragma once

#include "ofMain.h"
#include "ofxMol.h"

class ofApp : public ofBaseApp
{

public:
    void setup();
    void draw();
    void update();
    
protected:
    OfxMol::System system;
    ofMesh mesh;
};
