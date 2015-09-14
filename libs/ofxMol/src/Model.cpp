// Copyright (c) 2015 Davide Rambaldi.
// All rights reserved.
//
// This file is part of ofxMol.
//
// ofxMol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ofxMol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ofxMol.  If not, see <http://www.gnu.org/licenses/>.
//
//
// Additional permission under GNU GPL version 3 section 7
//
// If you modify this Library, or any covered work, by linking or
// combining it with ofxMol (or a modified version of that library), the
// licensors of this Library grant you additional permission to convey
// the resulting work. Corresponding Source for a non-source form of
// such a combination shall include the source code for the parts of CGAL
// used as well as that of the covered work.
//
//
//
// Author(s)     :  Davide Rambaldi

#include "ofxMol/Model.h"

namespace OfxMol
{
    Model::Model(): _model_number(0)
    {
        atoms.clear();
        coarse_atoms.clear();
    }
    
    Model::Model(int nm) : _model_number(nm)
    {
        atoms.clear();
        coarse_atoms.clear();
    }
    
    Model::~Model()
    {
        atoms.clear();
        coarse_atoms.clear();
    }
    
    void Model::add_atom(OfxMol::Atom &atom)
    {
        atoms.push_back(atom);
    }
    
    void Model::add_coarse_atom(OfxMol::Coarse_Atom &atom)
    {
        coarse_atoms.push_back(atom);
    }
    
    std::string Model::log()
    {
        std::string msg = "Model number: " + ofToString( model_number()) + "\n";
        msg += "Atoms: " + ofToString(number_of_atoms()) + "\n";
        msg += "Coarse atoms: " + ofToString(number_of_coarse_atoms()) + "\n";
        return msg;
    }
    
    //! Create a polyline on backbone
    ofPolyline Model::backbonePoly()
    {
        ofPolyline line;
        
        for (Const_atoms_iterator atm=atoms_begin(); atm!=atoms_end(); ++atm)
        {
            if (atm->is_backbone())
            {
                line.addVertex(atm->position());
            }
        }
        
        return line;
    }
    
    ofMesh Model::coarseAtomsMesh(int resolution)
    {
        ofMesh mesh;
        mesh.clear();
        mesh.enableColors();
        
        for (Coarse_atoms_iterator atm=coarse_atoms_begin(); atm!=coarse_atoms_end(); ++atm)
        {
            ofSpherePrimitive sphere = atm->sphere(resolution);
            vector<ofMeshFace> triangles = sphere.getMesh().getUniqueFaces();
            updateMesh(mesh, triangles, atm->position(),atm->getColor());
        }

        return mesh;
    }
    
    ofMesh Model::coarseAtomsMesh(ofColor color, int resolution)
    {
        ofMesh mesh;
        mesh.clear();
        mesh.enableColors();
        
        for (Coarse_atoms_iterator atm=coarse_atoms_begin(); atm!=coarse_atoms_end(); ++atm)
        {
            ofSpherePrimitive sphere = atm->sphere(resolution);
            vector<ofMeshFace> triangles = sphere.getMesh().getUniqueFaces();
            updateMesh(mesh, triangles, atm->position(), color);
        }
        
        return mesh;
    }
    
    //! Create a sphere for each atom in the system.
    //! Radius: given that default radius id 1.8 (check ESBTL code)
    //! I have 2 methods: atomsMesh with no radius arg use atom radius.
    //! The method with radius arg override arom radius.
    ofMesh Model::atomsMesh(int resolution)
    {
        ofMesh mesh;
        mesh.clear();
        
        for (Atoms_iterator atm=atoms_begin(); atm!=atoms_end(); ++atm)
        {
            ofSpherePrimitive sphere = atm->sphere(resolution);
            vector<ofMeshFace> triangles = sphere.getMesh().getUniqueFaces();
            updateMesh(mesh, triangles, atm->position());
        }
        return mesh;
    }
    
    // ! method for mesh gen with an arbitrary radius for all atoms
    ofMesh Model::atomsMesh(float radius, int resolution)
    {
        ofMesh mesh;
        mesh.clear();
        
        for (Atoms_iterator atm=atoms_begin(); atm!=atoms_end(); ++atm)
        {
            ofSpherePrimitive sphere = atm->sphere(radius, resolution);
            vector<ofMeshFace> triangles = sphere.getMesh().getUniqueFaces();
            updateMesh(mesh, triangles, atm->position());
        }
        return mesh;
    }
    
    void Model::updateMesh(ofMesh& mesh, const vector<ofMeshFace> &triangles,  const ofVec3f position)
    {
        for (int j=0; j<triangles.size(); j++)
        {
            ofMeshFace face = triangles[j];
            mesh.addVertex(face.getVertex(0) + position);
            mesh.addNormal(face.getNormal(0));
            mesh.addVertex(face.getVertex(1) + position);
            mesh.addNormal(face.getNormal(1));
            mesh.addVertex(face.getVertex(2) + position);
            mesh.addNormal(face.getNormal(2));
        }
    }
    
    void Model::updateMesh(ofMesh& mesh, const vector<ofMeshFace> &triangles,  const ofVec3f position, const ofColor color)
    {
        for (int j=0; j<triangles.size(); j++)
        {
            ofMeshFace face = triangles[j];
            
            mesh.addVertex(face.getVertex(0) + position);
            mesh.addNormal(face.getNormal(0));
            mesh.addColor(color);
            
            mesh.addVertex(face.getVertex(1) + position);
            mesh.addNormal(face.getNormal(1));
            mesh.addColor(color);
            
            mesh.addVertex(face.getVertex(2) + position);
            mesh.addNormal(face.getNormal(2));
            mesh.addColor(color);
        }
    }

    //! create a point cloud using atoms and atoms colors
    ofMesh Model::atomsPointCloud()
    {
        ofMesh mesh;
        mesh.enableColors();
        
        for (Const_atoms_iterator atm=atoms_begin(); atm!=atoms_end(); ++atm)
        {
            mesh.addColor(atm->getColor());
            mesh.addVertex(atm->position());
        }
        
        return mesh;
    }
}





