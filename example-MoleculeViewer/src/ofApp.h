#pragma once

#include "ofMain.h"
#include "ofxMol.h"

#define DATADIR "pdb"

class ofApp : public ofBaseApp
{

public:
    void setup();
    void update();
    void draw();

    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    void drawInteractionArea();
    
    void loadMolecule(ofFile file);
    
protected:
    bool bShowHelp;
    bool bCoarseAtomsMesh;
    bool bAtomsPointCloud;
    bool bAtoms;
    bool bBackbone;
    
    ofEasyCam cam;
    ofLight pointLightLeft;
    ofLight pointLightRight;
    ofMaterial material;

    OfxMol::System system;
    ofMesh cAtmMesh;
    ofMesh atomsCloud;
    ofMesh atoms;
    ofPolyline backbone;
    
    ofDirectory dataDir;
    std::vector<ofFile> pdbFiles;
    int currentFile;
};
