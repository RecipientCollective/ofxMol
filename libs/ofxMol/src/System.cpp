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

#include "ofxMol/System.h"

namespace OfxMol
{
    // SYSTEM IS DEFAULT WITH COARSE GRAIN
    typedef ESBTL::All_atom_system_builder<ESBTL::Default_system_with_coarse_grain> OfxMol_Builder;
    typedef ESBTL::Coarse_atoms_iterators<ESBTL::Default_system_with_coarse_grain::Model>::iterator OfxMol_Coarse_atoms_iterator;
    
    // OCCUPANCY POLICIES
    typedef ESBTL::Accept_all_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_all_occupancy_policy;
    typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;
    
    void System::setupSimple(std::string &path)
    {
        // simple line selector: all atoms and hetero-atoms are in the one system.
        ESBTL::PDB_line_selector sel;
        /*
         * Class responsible for building a coarse grain model
         * Creator of a coarse grain model creating up to two pseudo-atoms:
         * one as the barycenter of backbone atoms (if relevant) and one as the barycenter
         * of non-barycenter atoms.
         * Backbone atoms are identified using the global function ESBTL::is_backbone.
         */
        ESBTL::Coarse_creator_two_barycenters<ESBTL::Default_system_with_coarse_grain::Residue> creator;
        
        systems.clear();
        models.clear();
        water_models.clear();
        
        //Build the system from the PDB file.
        OfxMol_Builder builder(systems,sel.max_nb_systems());
        
        if (ESBTL::read_a_pdb_file(path,sel,builder,Accept_all_occupancy_policy()))
        {
            if (systems.empty() || systems.size() != 1)
            {
                ofLogFatalError() << "[ofxMol::System] Can create systems from input file: " << path << " using Mode SIMPLE";
            }
            
            if (systems[0].has_no_model())
            {
                ofLogError() << "[ofxMol::System] No models found in file: " << path;
            }
            
            // BEGIN MODEL LOOP
            for (ESBTL::Default_system_with_coarse_grain::Models_iterator it_model=systems[0].models_begin(); it_model!=systems[0].models_end(); ++it_model)
            {
                // create coarse atoms in models
                ESBTL::Default_system_with_coarse_grain::Model& model =* it_model;
                for (ESBTL::Default_system_with_coarse_grain::Model::Residues_iterator it_res=model.residues_begin(); it_res!=model.residues_end(); ++it_res)
                {
                    it_res->create_coarse_atoms(creator);
                }
                
                // ofxModel
                OfxMol::Model ofxModel = OfxMol::Model(model.model_number());
                
                for (ESBTL::Default_system_with_coarse_grain::Model::Atoms_iterator it_atm=model.atoms_begin(); it_atm!=model.atoms_end(); ++it_atm)
                {
                    // add atom
                    OfxMol::Atom atom = OfxMol::Atom(*it_atm);
                    // ofLogNotice() << it_atm->temperature_factor();
                    ofxModel.add_atom(atom);
                }
                
                for (OfxMol_Coarse_atoms_iterator itc=ESBTL::coarse_atoms_begin(model); itc!=ESBTL::coarse_atoms_end(model); ++itc)
                {
                    // add coarse atoms
                    OfxMol::Coarse_Atom coarse_atom = OfxMol::Coarse_Atom(*itc);
                    ofxModel.add_coarse_atom(coarse_atom);
                }
                
                // push model
                models.push_back(ofxModel);
            }
            // END MODEL LOOP
            
            ofLogNotice() << "[ofxMol::System] Setup complete for file: " << path;
        }
        else
        {
            ofLogFatalError() << "[ofxMol::System] Setup incomplete for file: " << path;
        }
    }
    
    void System::setupAdvanced(std::string &path)
    {
        /** This is a line selector defining two systems:
         * - the first one contains all heavy atoms that do not belong to a water molecule.
         * - the second one contains all containing heavy atoms of water molecules.
         */
        ESBTL::PDB_line_selector_two_systems sel;
        /*
         * Class responsible for building a coarse grain model
         * Creator of a coarse grain model creating up to two pseudo-atoms:
         * one as the barycenter of backbone atoms (if relevant) and one as the barycenter
         * of non-barycenter atoms.
         * Backbone atoms are identified using the global function ESBTL::is_backbone.
         */
        ESBTL::Coarse_creator_two_barycenters<ESBTL::Default_system_with_coarse_grain::Residue> creator;
        
        systems.clear();
        models.clear();
        water_models.clear();
        
        //Build the system from the PDB file.
        OfxMol_Builder builder(systems,sel.max_nb_systems());
        
        if (ESBTL::read_a_pdb_file(path,sel,builder,Accept_all_occupancy_policy()))
        {
            if (systems.empty() || systems.size() != 2)
            {
                ofLogFatalError() << "[ofxMol::System] Can create systems from input file: " << path << " using Mode ADVANCED";
            }
            
            if (systems[0].has_no_model())
            {
                ofLogError() << "[ofxMol::System] No atoms found in file: " << path;
            }
            
            systems[0].name() = "atoms";
            systems[1].name() = "water";
            
            // BEGIN MODEL LOOP
            for (ESBTL::Default_system_with_coarse_grain::Models_iterator it_model=systems[0].models_begin(); it_model!=systems[0].models_end(); ++it_model)
            {
                // create coarse atoms in models
                ESBTL::Default_system_with_coarse_grain::Model& model =* it_model;
                for (ESBTL::Default_system_with_coarse_grain::Model::Residues_iterator it_res=model.residues_begin(); it_res!=model.residues_end(); ++it_res)
                {
                    it_res->create_coarse_atoms(creator);
                }
                
                // ofxModel
                OfxMol::Model ofxModel = OfxMol::Model(model.model_number());
                
                for (ESBTL::Default_system_with_coarse_grain::Model::Atoms_iterator it_atm=model.atoms_begin(); it_atm!=model.atoms_end(); ++it_atm)
                {
                    // add atom
                    OfxMol::Atom atom = OfxMol::Atom(*it_atm);
                    // ofLogNotice() << it_atm->temperature_factor();
                    ofxModel.add_atom(atom);
                }
                
                for (OfxMol_Coarse_atoms_iterator itc=ESBTL::coarse_atoms_begin(model); itc!=ESBTL::coarse_atoms_end(model); ++itc)
                {
                    // add coarse atoms
                    OfxMol::Coarse_Atom coarse_atom = OfxMol::Coarse_Atom(*itc);
                    ofxModel.add_coarse_atom(coarse_atom);
                }
                
                // push model
                models.push_back(ofxModel);
            }
            // END MODEL LOOP
            
            
            // BEGIN WATER MODEL LOOP
            for (ESBTL::Default_system_with_coarse_grain::Models_iterator it_model=systems[1].models_begin(); it_model!=systems[1].models_end(); ++it_model)
            {
                const ESBTL::Default_system_with_coarse_grain::Model& model =* it_model;
                OfxMol::Model ofxWaterModel = OfxMol::Model(model.model_number());
                
                for (ESBTL::Default_system_with_coarse_grain::Model::Atoms_const_iterator it_atm=model.atoms_begin(); it_atm!=model.atoms_end(); ++it_atm)
                {
                    OfxMol::Atom atom = OfxMol::Atom(*it_atm);
                    ofxWaterModel.add_atom(atom);
                }
                
                water_models.push_back(ofxWaterModel);
            }
            // BEGIN WATER MODEL LOOP
            ofLogNotice() << "[ofxMol::System] Setup complete for file: " << path;
        }
        else
        {
            ofLogFatalError() << "[ofxMol::System] Setup incomplete for file: " << path;
        }
    }
    
    void System::setup(std::string path, SetupMode mode)
    {
        switch (mode)
        {
            case SIMPLE:
                setupSimple(path);
                break;
            case ADVANCED:
                setupAdvanced(path);
                break;
            default:
                break;
        }

    }
}