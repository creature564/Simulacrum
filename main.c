//
//  sidr2.c
//  SIDR v2
//
//  Created by Adam Case on 01/31/2021.
//  Copyright © 2021 Adam Case.
//  All rights reserved.
//

#include <stdio.h>

#include "pipeline.h"
#include "param.h"
#include "pylink.h"


int main(int argc, char *argv[]) {

    PARAM param;
    initPARAM(&param, argc, argv);

    SEQDATA *data = run_pipeline(&param);

    displaySEQDATA(data);

    printf("analysis: %d\n", run_analysis(&param, data));

    freeSEQDATA(data);
    freePARAM(&param);

    return 0;

}
