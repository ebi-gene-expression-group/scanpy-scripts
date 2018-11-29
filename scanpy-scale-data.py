#!/usr/bin/env python

from __future__ import print_function
import logging
import scanpy.api as sc
from scanpy_wrapper_utils import ScanpyArgParser, comma_separated_list
from scanpy_wrapper_utils import read_input_object, write_output_object


def main(args):
    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)

    if args.do_log:
        sc.pp.log1p(adata)

    if args.var_to_regress:
        sc.pp.regress_out(adata, args.var_to_regress)

    sc.pp.scale(adata, zero_center=args.zero_center, max_value=args.scale_max)

    write_output_object(adata, args.output_object_file, args.output_format)

    print(adata)
    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Scale data to equal variance')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('-g', '--do-log',
                           action='store_true',
                           help='Log-transform the data: log(x+1)')
    argparser.add_argument('-V', '--var-to-regress',
                           type=comma_separated_list('var-to-regress', str),
                           default=[],
                           help='Variables to regress out')
    argparser.add_argument('-z', '--zero-center',
                           action='store_true',
                           dest='zero_center',
                           default=None,
                           help='Zero-center the data. '
                                'See also --no-zero-center. '
                                'If neither --zero-center nor --no-zero-center is set, '
                                'automatically set to FALSE if input is sparse, otherwise TRUE.')
    argparser.add_argument('-Z', '--no-zero-center',
                           action='store_false',
                           dest='zero_center',
                           help='Do not zero-center the data to handle sparse data efficiently. '
                                'See also --zero-center.')
    argparser.add_argument('-x', '--scale-max',
                           type=float,
                           default=None,
                           help='Truncate to this value after scaling. '
                                'Do not truncate if None. Default: None')
    args = argparser.get_args()

    main(args)
