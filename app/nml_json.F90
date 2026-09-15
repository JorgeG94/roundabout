!! Roundabout namelist knob-reference JSON emitter (P4 generator input).
program nml_json
   !! Host-only CLI tool: build a default `config_t`, register the full
   !! strict schema against it, and render every group + knob to JSON via
   !! `render_json`.  Verbatim sibling of `app/nml_doc.F90` -- it parses
   !! nothing, it introspects the live schema object.  Output path is
   !! arg 1 (default `nml_knobs.json`).
   !!
   !! `tools/gen_python_config.py` reads this JSON and emits the checked-in
   !! `python/rdb/_config_generated.py`.  `tests/test_python_config_drift`
   !! (python/tests/test_config_drift.py) regenerates both and diffs
   !! against the checked-in file -- this program is the first stage of
   !! that pipeline, so a drift-test regeneration command starts here:
   !!
   !!     ./rdb_nml_json tmp_local_artifacts/schema.json
   !!     python3 tools/gen_python_config.py tmp_local_artifacts/schema.json python/rdb/_config_generated.py
   use rdb_config, only: config_t, build_rdb_schema
   use rdb_nml_schema, only: nml_schema_t
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none

   type(config_t), target :: cfg
   type(nml_schema_t) :: schema
   character(len=512) :: out_path
   integer :: arg_len, io_stat, ngroups, nknobs, ig

   if (command_argument_count() >= 1) then
      call get_command_argument(1, out_path, arg_len, io_stat)
   else
      out_path = "nml_knobs.json"
   end if

   call build_rdb_schema(cfg, schema)
   call schema%render_json(trim(out_path))

   ngroups = size(schema%groups)
   nknobs = 0
   do ig = 1, ngroups
      if (allocated(schema%groups(ig)%keys)) nknobs = nknobs + size(schema%groups(ig)%keys)
   end do

   call logger%info("Wrote namelist knob JSON to "//trim(out_path)// &
                    " ("//to_char(ngroups)//" groups, "//to_char(nknobs)//" knobs)")

end program nml_json
