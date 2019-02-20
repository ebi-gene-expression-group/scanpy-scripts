#!/usr/bin/env python

from __future__ import print_function
import logging
import scanpy.api as sc
from scanpy_wrapper_utils import ScanpyArgParser, read_input_object, write_output_object
import scipy
import pandas as pd
import os
import shutil

def main(args):
    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)
    print(adata)

    if args.matrix_format == 'mtx_zip':

        # Assemble a sparse matrix
        mat = scipy.sparse.coo_matrix.transpose(scipy.sparse.coo_matrix(adata.X))
        n, m = mat.shape
        l = len(mat.data)
        header = '%%MatrixMarket matrix coordinate real general\n%\n{} {} {}\n'.format(n, m, l)
        df = pd.DataFrame({'row': 1+mat.row, 'col': 1+mat.col, 'data':mat.data}) 
        
        # Define the output directory from the zip file name
        outdir = os.path.splitext(args.output_matrix)[0]
        os.mkdir(outdir)

        # Write the component files
        f = open(outdir + '/matrix.mtx', 'a')
        f.write(header)
        df.to_csv(f, sep=' ', header=False, index=False)
        f.close()
       
        # Write row and column labels
        with open(outdir + '/barcodes.tsv', 'w') as f:
            for item in list(adata.obs.index):
                f.write("%s\n" % item)

        with open(outdir + '/genes.tsv', 'w') as f:
            for item in list(adata.var.index):
                f.write("%s\t%s\n" % (item, item) )

        # Create the archive
        shutil.make_archive(outdir, 'zip', '.', outdir)

        # Clean up
        shutil.rmtree(outdir)

        logging.info('Done')
        return 0

    else: 

        logging.error("Only 'mtx_zip' currently valid as output format")
        ys.exit(1)


if __name__ == '__main__':
    argparser = ScanpyArgParser()
    argparser.add_input_object()
    argparser.add_argument('-m', '--matrix_format',
                           choices=['mtx_zip'],
                           default='mtx_zip',
                           help='Format to store output matrix.')

    argparser.add_argument('-o', '--output_matrix',
                           default=None,
                           required=True,
                           help='File to store output matrix.')

    args = argparser.get_args()

    main(args)
