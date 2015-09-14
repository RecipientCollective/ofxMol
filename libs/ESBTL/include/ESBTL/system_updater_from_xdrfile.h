// Copyright (c) 2009-2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
//This file is part of ESBTL.
//
//ESBTL is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//ESBTL is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with ESBTL.  If not, see <http://www.gnu.org/licenses/>.
//
//
//Additional permission under GNU GPL version 3 section 7
//
//If you modify this Library, or any covered work, by linking or
//combining it with CGAL (or a modified version of that library), the
//licensors of this Library grant you additional permission to convey
//the resulting work. Corresponding Source for a non-source form of
//such a combination shall include the source code for the parts of CGAL
//used as well as that of the covered work. 
//
//
//
// Author(s)     :  Sébastien Loriot



#ifndef SYSTEM_UPDATER_FROM_XDRFILE_H
#define SYSTEM_UPDATER_FROM_XDRFILE_H

// Require library xdrfile available at ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz
// and refer to http://wiki.gromacs.org/index.php/XTC_Library
//For more info on how to read the xtc file format refer to
// http://www.gromacs.org/documentation/reference_4.0/online/xtc.html


#include <cstdlib>
#include <boost/tuple/tuple.hpp>
#include <iostream>

#define CPLUSPLUS
#include <xdrfile/xdrfile.h>
#undef CPLUSPLUS



namespace ESBTL{
  
  /**
    * Helper class responsible for updating the coordinates of a system
    * according to a trajectory read from an xtc file. The system must already
    * have been constructed (from a PDB file for example).
    *
    * Note that using this class requires that the xdrfile library and header to be 
    * installed on your system (available as a package for your distribution or 
    * at ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz).
    *
    * Refer to http://www.gromacs.org for more information.
    * \tparam System is the type of the system used.
    */
  template <class System>
  class System_updater_from_xdrfile{
    typedef std::map<unsigned,typename System::Atom*> Map_id_to_atom;
    Map_id_to_atom selected_atoms_;
    
    XDRFILE* input_file_;
    unsigned first_frame_id_;
    unsigned last_frame_id_;
    unsigned nb_atoms_to_read_;
    unsigned current_frame_;
  public:
    
    /**
      * Constructor. It requires an already build model (from a PDB file for example.)
      * \tparam System_iterator is an iterator over systems used.
      * \param begin is an iterator over the first system that need to be updated.
      * \param end is a past iterator over the last system that need to be updated.
      * \param model_selected_ indicates the number of the model to be updated.
      * \param xtc_fname is the path to the xtc file used to update the system coordinates.
      * \param max_atoms is an upper bound on the number of atoms contains in a frame.
      * \param read_from is the index of the first frame to be loaded.
      * \param read_to is the index of the last frame to be loaded.
      */
    template <class System_iterator>
    System_updater_from_xdrfile(System_iterator begin, System_iterator end,unsigned model_selected_,const std::string& xtc_fname,
                                unsigned max_atoms,unsigned read_from=1,
                                unsigned read_to=std::numeric_limits<unsigned>::max() ):first_frame_id_(read_from),last_frame_id_(read_to),nb_atoms_to_read_(max_atoms),current_frame_(0)
    {
      input_file_=xdrfile_open(xtc_fname.c_str(), "r");
      if (input_file_ == NULL){
        std::cerr << "Error while opening xtc file " << xtc_fname << std::endl;
        exit (EXIT_FAILURE);
      }
      
      for (System_iterator it_sys=begin;it_sys!=end;++it_sys){
        typename System::Model& model=it_sys->get_model(model_selected_);
        for (typename System::Model::Atoms_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
          assert( selected_atoms_.find(it_atm->atom_serial_number())==selected_atoms_.end() );
          selected_atoms_.insert(std::make_pair(it_atm->atom_serial_number(),&(*it_atm)));
        }
      }
      
      //check the file contains at least one frame
      int magic=0;
      int result=xdrfile_read_int(&magic,1,input_file_);
      if (result==0){
        std::cerr << "xtc file provided contains no frame: error while reading magic number" << std::endl;
        exit (EXIT_FAILURE);
      }
      assert (magic==1995); //checks XDR magic number    
      
      if (first_frame_id_!=1)
        next_frame(true);
      
      assert(current_frame_+1==first_frame_id_);
    }
    
    ~System_updater_from_xdrfile(){
      //only when reading stops at a given frame id
      if (input_file_!=NULL)
        xdrfile_close(input_file_);
    }
      
    
    //return true if next_frame can be loaded
    //return in addition the corresponding simulation step and the simulation time.
    /**
      * Loads the next frame into the system (update coordinates).
      * \return the tuple returned contains a boolean (indicating whether the next frame could be loaded),
      * a double (indicating the time of the frame in the simulation) and an integer (indicating the simulation step).
      */
    boost::tuple<bool,double,int> next_frame(bool is_init=false){
      if (!has_more_frames())
        return boost::make_tuple(false,-1,-1);
      
      int magic,result, natoms,step;
      float time, prec, box[9];//, min[3], max[3];
      float xyz[3*nb_atoms_to_read_];

      do
      {
        ++current_frame_;
       // Read the compressed xtc file
        xdrfile_read_int(&natoms,1,input_file_); // reads number of atoms in the frame
        assert(natoms == static_cast<int>(nb_atoms_to_read_)); // number of atoms in the PDB record matches the frame
        xdrfile_read_int(&step,1,input_file_); // reads the step number
        xdrfile_read_float(&time,1,input_file_);
        xdrfile_read_float(box,9,input_file_); // reads the box dimensions
        result = xdrfile_decompress_coord_float(xyz,&natoms,&prec,input_file_); //reads the coordinates in xyz
        assert(result != 0);
        
        result=xdrfile_read_int(&magic,1,input_file_);
        if (result==0){
          if (current_frame_ <= first_frame_id_ - 1){
            std::cerr << "Starting frame id is greater than the number of frames in the provided file" << std::endl;
            exit(EXIT_FAILURE);
          }
          xdrfile_close(input_file_);
          input_file_=NULL;
        }
        assert (result==0 || magic==1995); //checks XDR magic number
      }
      while( current_frame_ < first_frame_id_ - 1 );
      
      //update the coordinates
      if (!is_init)
        for(typename Map_id_to_atom::iterator it_sel= selected_atoms_.begin(); it_sel!=selected_atoms_.end();++it_sel){
          unsigned i=it_sel->first;
          typename System::Atom::Point_3 new_center(
            xyz[3*(i-1)]*10, 
            xyz[3*(i-1)+1]*10,
            xyz[3*(i-1)+2]*10
          );//*10 because gromacs is in nanometer

          static_cast<typename System::Atom::Point_3&>(*it_sel->second)=new_center;
        }
      return boost::make_tuple(true,time,step);
    }
    
    /** Indicates whether a frame is available to update the system.*/
    bool has_more_frames(){
      return (input_file_!=NULL && current_frame_ !=last_frame_id_);
    }
    
    const unsigned& current_frame_id(){return current_frame_;}
  };
  
}//namespace ESBTL


















#endif //SYSTEM_UPDATER_FROM_XDRFILE_H

