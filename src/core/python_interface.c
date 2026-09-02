#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Declare the existing C main function
extern int main(int argc, char *argv[]);

static PyObject* py_wrapper(PyObject *self, PyObject *args) {
    PyObject *list_obj;

    // Expect a single list argument from Python
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &list_obj)) {
        return NULL;
    }

    Py_ssize_t c_argc = PyList_Size(list_obj);
    char **c_args = (char **)malloc(c_argc * sizeof(char *));
    if (!c_args) {
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < c_argc; i++) {
        PyObject *item = PyList_GetItem(list_obj, i);
        if (!PyUnicode_Check(item)) {
            // Cleanup and raise error if not string
            free(c_args);
            PyErr_SetString(PyExc_TypeError, "All arguments must be strings.");
            return NULL;
        }
        c_args[i] = (char *)PyUnicode_AsUTF8(item);
    }

    // Call the original C main function directly with the parsed arguments!
    int result = main((int)c_argc, c_args);

    free(c_args);

    if (result != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Simulation execution failed inside the C core.");
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyMethodDef CoreMethods[] = {
    {"wrapper", py_wrapper, METH_VARARGS, "Execute the simulation engine with arguments."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef coremodule = {
    PyModuleDef_HEAD_INIT,
    "_core",
    NULL,
    -1,
    CoreMethods
};

PyMODINIT_FUNC PyInit__core(void) {
    return PyModule_Create(&coremodule);
}