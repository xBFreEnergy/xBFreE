# ##############################################################################
#                            GPLv3 LICENSE INFO                                #
#  Copyright (c) 2023.  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco  #
#                                                                              #
#  This program is free software; you can redistribute it and/or modify it     #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

import sys
from os.path import split
import logging

from xBFreE.mmpbsa.main import MMPBSA_App
from xBFreE.exceptions import InputError
from xBFreE.mmpbsa.infofile import InfoFile

def mmpbsa(parser, MPI):
    # Instantiate the main MMPBSA_App
    app = MMPBSA_App(MPI)
    # Read the command-line arguments
    app.get_parser(parser)

    # If we're not rewriting output do whole shebang, otherwise load info and parms
    # Throw up a barrier before and after running the actual calcs
    if not app.FILES.rewrite_output:
        try:
            app.read_input_file()
        except InputError as e:
            sys.stderr.write('%s: %s' % (type(e).__name__, e) + '\n')
            sys.stderr.write('  Enter `%s --help` for help\n' %
                             (split(sys.argv[0])[1]))
            sys.exit(1)
        app.process_input()
        app.check_for_bad_input()
        app.make_prmtops()
        app.loadcheck_prmtops()
        app.file_setup()
        app.run_mmpbsa()
    # If we are rewriting output, load the info and check prmtops
    else:
        info = InfoFile(app, True)
        info.read_info()
        app.loadcheck_prmtops()

    # Now we parse the output, print, and finish
    app.parse_output_files()
    app.write_final_outputs()
    app.finalize()


if __name__ == '__main__':


    logging.info('Finished')

    # mmpbsa()
