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

#include <ESBTL/default.h>

#include "ofMain.h"
#include "ofxMol/Model.h"

namespace OfxMol
{
    // modes:
    // SIMPLE: all atoms in one system
    //
    // ADVANCED:
    //  2 systems the first one contains all heavy atoms that do not belong to a water molecule.
    //  the second one contains all containing heavy atoms of water molecules. Hydrogen atoms are discarded
    //  Use water_models_begin and water_models_end to iterate over water models in second system
    //
    // CHAINS: by chains (todo)
    //
    enum SetupMode
    {
        SIMPLE,
        ADVANCED
    };
    
    class System
    {
    public:
        void setup(std::string path, SetupMode mode = SIMPLE);
        
        //! iterator for models
        typedef std::vector<OfxMol::Model>::const_iterator Const_models_iterator;
        typedef std::vector<OfxMol::Model>::iterator Models_iterator;
        
        Const_models_iterator models_begin() const
        {
            return models.begin();
        }
        
        Const_models_iterator models_end() const
        {
            return models.end();
        }
        
        Models_iterator models_begin()
        {
            return models.begin();
        }
        
        Models_iterator models_end()
        {
            return models.end();
        }
        
        Model &getModel(unsigned int i)
        {
            return models[i];
        }
        
        unsigned int number_of_models() const
        {
            return models.size();
        }
        
        
        //! iterator for water models
        Const_models_iterator water_models_begin() const
        {
            return water_models.begin();
        }
        
        Const_models_iterator water_models_end() const
        {
            return water_models.end();
        }
        
        Models_iterator water_models_begin()
        {
            return water_models.begin();
        }
        
        Models_iterator water_models_end()
        {
            return water_models.end();
        }
        
        Model getWaterModel(unsigned int i)
        {
            return water_models[i];
        }
        
        unsigned int number_of_water_models() const
        {
            return water_models.size();
        }
        
    protected:
        void setupSimple(std::string &path);
        void setupAdvanced(std::string &path);
        std::vector<ESBTL::Default_system_with_coarse_grain> systems;
        std::vector<OfxMol::Model> models;
        std::vector<OfxMol::Model> water_models;
    };
}