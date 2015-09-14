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

#pragma once

#include "ofxMol/Atom.h"
#include "ofxMol/Coarse_Atom.h"

namespace OfxMol
{
    class Model
    {
    public:
        
        Model();
        Model(int nm);
        ~Model();
        
        inline const int model_number(){ return _model_number; }
        
        //! Add atoms
        void add_atom(OfxMol::Atom &atom);
        inline const int number_of_atoms() { return atoms.size(); }
        void add_coarse_atom(OfxMol::Coarse_Atom &atom);
        inline const int number_of_coarse_atoms() { return coarse_atoms.size(); }
        
        //! Logger
        std::string log();
        
        //! iterators for atoms
        typedef std::vector<OfxMol::Atom>::const_iterator Const_atoms_iterator;
        typedef std::vector<OfxMol::Atom>::iterator Atoms_iterator;
        
        inline Const_atoms_iterator atoms_begin() const
        {
            return atoms.begin();
        }
        
        inline Const_atoms_iterator atoms_end() const
        {
            return atoms.end();
        }
        
        inline Atoms_iterator atoms_begin()
        {
            return atoms.begin();
        }
        
        inline Atoms_iterator atoms_end()
        {
            return atoms.end();
        }
        
        Atom getAtom(unsigned int i)
        {
            return atoms[i];
        }
        
        //! iterators for coarse atoms
        typedef std::vector<OfxMol::Coarse_Atom>::const_iterator Const_coarse_atoms_iterator;
        typedef std::vector<OfxMol::Coarse_Atom>::iterator Coarse_atoms_iterator;
        
        inline Const_coarse_atoms_iterator coarse_atoms_begin() const
        {
            return coarse_atoms.begin();
        }
        
        inline Const_coarse_atoms_iterator coarse_atoms_end() const
        {
            return coarse_atoms.end();
        }
        
        inline Coarse_atoms_iterator coarse_atoms_begin()
        {
            return coarse_atoms.begin();
        }
        
        inline Coarse_atoms_iterator coarse_atoms_end()
        {
            return coarse_atoms.end();
        }
        
        Coarse_Atom getCoarseAtom(unsigned int i)
        {
            return coarse_atoms[i];
        }
        
        //! generators
        ofMesh atomsPointCloud();
        ofMesh atomsMesh(int resolution = 16);
        ofMesh atomsMesh(float radius, int resolution = 16);
        ofMesh coarseAtomsMesh(ofColor color, int resolution = 16);
        ofMesh coarseAtomsMesh(int resolution = 16);
        ofPolyline backbonePoly();
        
        
    protected:
        int _model_number;
        std::vector<OfxMol::Atom> atoms; // atoms
        std::vector<OfxMol::Coarse_Atom> coarse_atoms; // coarse atoms
        void updateMesh(ofMesh& mesh, const vector<ofMeshFace> &triangles, const ofVec3f position, const ofColor color);
        
         void updateMesh(ofMesh& mesh, const vector<ofMeshFace> &triangles, const ofVec3f position);
    };
}





