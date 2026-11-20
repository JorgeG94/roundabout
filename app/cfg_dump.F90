!! Resolved-configuration dump: the equivalence oracle between input formats.
program cfg_dump
   !! `rdb_cfg_dump <input> <out>`: read `<input>` (namelist, or TOML for a
   !! `.toml` name) through `read_config`, then write the resolved
   !! configuration -- every group and key, post-cascade -- to `<out>`
   !! (`write_doc_all`) and, at full real precision, to `<out>.exact`
   !! (`write_resolved`; the doc form rounds reals to 11 digits).
   !!
   !! Both files are a function of the resolved configuration only: no path,
   !! timestamp or input-format text reaches them, so an `.nml` and the
   !! `.toml` that spells the same configuration give byte-identical output.
   !! Exit 0 on success, 1 when the input is missing or refused.
   use rdb_config, only: config_t, read_config
   use rdb_nml_schema, only: nml_schema_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use pic_logger, only: logger => global_logger
   implicit none

   type(config_t), target :: cfg
   type(nml_schema_t) :: schema
   character(len=1024) :: in_path, out_path
   integer :: ierr
   logical :: exists

   if (command_argument_count() < 2) then
      call logger%error("usage: rdb_cfg_dump <input.nml|input.toml> <out>")
      stop 2
   end if
   call get_command_argument(1, in_path)
   call get_command_argument(2, out_path)

   inquire (file=trim(in_path), exist=exists)
   if (.not. exists) then
      call logger%error("rdb_cfg_dump: cannot open "//trim(in_path))
      stop 1
   end if

   call read_config(trim(in_path), cfg, schema, ierr)
   if (ierr /= OCEAN_STATUS_OK) then
      call logger%error("rdb_cfg_dump: configuration refused (see log)")
      stop 1
   end if
   call schema%write_doc_all(trim(out_path))
   call schema%write_resolved(trim(out_path)//".exact")
end program cfg_dump
