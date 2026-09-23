
# ── EasleeConfig — standalone run-configuration layer ─────────────────────────
# Loadable WITHOUT the heavy EASLEE/CUDA stack, so the GUI form, validation,
# CLI help, and FLEXPART import are fast.  build_run / run_spec! (which need
# EASLEE proper) live in easlee/build_run.jl and are included in Main by the
# scripts/CLI after runtime.jl — the runtime path is always
# read_runspec(TOML) → build_run / run_spec!.

module EasleeConfig

include("units.jl");      using .Units
include("runspec.jl");    using .RunSpecs
include("species_db.jl"); using .SpeciesDB
include("toml_io.jl");    using .ConfigTOML
include("schema.jl");     using .Schema
include("flexpart_import.jl"); using .FlexpartImport

# re-export the public surface
using .RunSpecs:  RunSpec, MetSource, Species, ReleaseBox, Physics, OutputSpec, validate_runspec
using .ConfigTOML: read_runspec, write_runspec, runspec_to_dict, dict_to_runspec
using .Schema:    schema_fields, schema_json
using .FlexpartImport: import_flexpart
using .Units:     parse_duration_s, parse_length_m, parse_mass_kg

export RunSpec, MetSource, Species, ReleaseBox, Physics, OutputSpec, validate_runspec
export read_runspec, write_runspec, runspec_to_dict, dict_to_runspec
export schema_fields, schema_json, import_flexpart

end # module EasleeConfig
