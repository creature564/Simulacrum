

// https://docs.python.org/3/extending/embedding.html
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "pylink.h"

int run_analysis(PARAM *param, SEQDATA *data) {

    Py_Initialize();
    PyRun_SimpleString("import sys\n"
                       "sys.path.insert(0, \"/jlf/acmccormack1/SIDR2.0/src/\")\n");
    PyObject *file = PyUnicode_DecodeFSDefault("analysis");
    PyObject *module = PyImport_Import(file);
    Py_DECREF(file);

    if (module != NULL) {

        PyObject *function = PyObject_GetAttrString(module, "run_analysis");

        if (function && PyCallable_Check(function)) {

            uint32_t regions = numRegions(data);
            uint32_t kmers = numKmers(data);
            PyObject *args = PyTuple_New(8);

            PyObject *seq_list = PyList_New(regions);
            PyObject *len_list = PyList_New(regions);
            PyObject *gc_list = PyList_New(regions);
            PyObject *cov_list = PyList_New(regions);
            PyObject *blast_list = PyList_New(regions);
            PyObject *tax_list = PyList_New(regions);
            PyObject *histFile_list = PyList_New(regions);
            PyObject *kmer_list = PyList_New(kmers);

            if (!seq_list || !len_list || !gc_list || !cov_list || !blast_list || !tax_list || !histFile_list || !kmer_list) {

                printf("Cannot create data lists\n");
                Py_FinalizeEx();
                return -1;

            }

            for(int i = 0; i < regions; i++) {

                PyList_SetItem(seq_list, i, Py_BuildValue("s", get_region_name(data, i)));
                PyList_SetItem(len_list, i, Py_BuildValue("I", get_region_length(data, i)));
                PyList_SetItem(gc_list, i, Py_BuildValue("f", get_region_gc(data, i)));
                PyList_SetItem(cov_list, i, Py_BuildValue("f", get_region_cov(data, i)));
                PyList_SetItem(blast_list, i, Py_BuildValue("i", get_region_blast(data, i)));
                PyList_SetItem(tax_list, i, Py_BuildValue("i", get_region_tax(data, i)));

                PyObject *regionHistFiles = PyList_New(kmers);
                for(int k = 0; k < kmers; k++) {
                    PyList_SetItem(regionHistFiles, k, Py_BuildValue("s", get_region_histfile(data, i, k)));
                }

                PyList_SetItem(histFile_list, i, regionHistFiles);

            }

            for(int k = 0; k < kmers; k++) {
                PyList_SetItem(kmer_list, k, Py_BuildValue("i", param->kmer_list[k]));
            }


            PyTuple_SetItem(args, 0, seq_list);
            PyTuple_SetItem(args, 1, len_list);
            PyTuple_SetItem(args, 2, gc_list);
            PyTuple_SetItem(args, 3, cov_list);
            PyTuple_SetItem(args, 4, blast_list);
            PyTuple_SetItem(args, 5, tax_list);
            PyTuple_SetItem(args, 6, histFile_list);
            PyTuple_SetItem(args, 7, kmer_list);


            PyObject *ret = PyObject_CallObject(function, args);
            Py_DECREF(args);

            if (ret != NULL) printf("Analysis successfully executed\n");

            else {

                PyErr_Print();
                printf("Call failed\n");

            }
        }

        else {

            if (PyErr_Occurred()) PyErr_Print();
            printf("Cannot find function\n");

        }

    }

    else {

        PyErr_Print();
        printf("Failed to load\n");

    }

    return Py_FinalizeEx();

}
