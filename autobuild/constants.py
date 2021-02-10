class Constants:
    EXECUTABLE = (
        ""
    )
    JOB = (
        ""
    )

    RHOFIT = (
        "#!/bin/bash \n"
        "conda activate pandda \n"
        "source ~/.bashrc \n"
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "{pandda_rhofit} -map {event_map} -mtz {mtz} -pdb {pdb} -cif {cif} -out {out_dir_path}"
    )

