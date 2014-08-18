#ifdef __cplusplus
extern "C" void zpptri_(
	const char &uplo,		// (input)
	const int &n,			// (input)
	complex<double> *ap,		// ap[n*(n+1)/2] (input/output)
	int &info			// (output)
	);
#else /* ! __cplusplus */
void zpptri_(
	const char *uplo,		/* (input) */
	const int *n,			/* (input) */
	complex_double *ap,		/* ap[n*(n+1)/2] (input/output) */
	int *info			/* (output) */
	);
#endif /* ! __cplusplus */

