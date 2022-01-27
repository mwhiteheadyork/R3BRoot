/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BPRPCREADER_H
#define R3BPRPCREADER_H

/* Include the parent reader class definition */
#include "R3BReader.h"

/* The structure containing the data layout of PTOF detector */
struct EXT_STR_h101_RPC_t;
typedef struct EXT_STR_h101_RPC_t EXT_STR_h101_RPC;
typedef struct EXT_STR_h101_RPC_onion_t EXT_STR_h101_RPC_onion;
class ext_data_struct_info;


/* Forward declaration of FairLogger */
class FairLogger;

/*
 * Reader class for PTOF detector
 */
class R3BRpcReader : public R3BReader
{
  public:
    /* Default constructor / destructor */
    R3BRpcReader(EXT_STR_h101_RPC*);
    ~R3BRpcReader();

    /*
     * Functions needing implementation:
     * These are pure virtual in the base class
     */

    /*
     * Initialise the detector specific data items
     * in the ext_data_struct_info list of items.
     * This is called by the R3BUcesbSource.
     */
    Bool_t Init(ext_data_struct_info*);

    /*
     * Copy data from the ucesb data stream to
     * native data containers, called by R3BUcesbSource.
     */
    Bool_t Read();

    /*
     * Reset the internal output array.
     */
    void Reset();

  private:
    /* Reader specific data structure from ucesb */
    EXT_STR_h101_RPC* fData;

    /*
     * This data offset describes where in the main
     * data structure the detector specific data for this
     * PTOF instance begins. Use e.g.:
     *   offsetof(EXT_STR_h101, PTOF);
     */
    UInt_t fOffset;

    /* FairLogger */
    FairLogger* fLogger;

    // An event counter
    unsigned int fNEvent;

    /*
     * An array where the read R3BPaddleTamexMappedData
     * items are stored.
     */
    TClonesArray* fArrayRef; /**< Output array. */
    TClonesArray* fArrayPmt; /**< Output array. */
    TClonesArray* fArrayStrip; /**< Output array. */

    Int_t fCoarseReference; // earliest coarse counter to detect coarse counter wrap

    /*
     * Helper functions to structure the reader code.
     */

  public:
    ClassDef(R3BRpcReader, 0);
};

#endif
