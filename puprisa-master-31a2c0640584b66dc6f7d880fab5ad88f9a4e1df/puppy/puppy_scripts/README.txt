Puppy, the puprisa batch processing engine.

All puprisa modules should have a batch option, and the ability to be fed an image stack, e.g.:

[euFraction, pheoFraction] = puprisa_linearSpectralDecomp( imageStack, refSpecFile, other options, 'batch' )

If the last argument is 'batch', the processing will execute immediately and return a result.