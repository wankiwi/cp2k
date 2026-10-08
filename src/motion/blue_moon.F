!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE blue_moon
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE blue_moon_types,                 ONLY: blue_moon_entry_type,&
                                              blue_moon_state_type
   USE cell_types,                      ONLY: cell_type
   USE colvar_methods,                  ONLY: colvar_eval_mol_f
   USE colvar_types,                    ONLY: &
        angle_colvar_id, coord_colvar_id, dist_colvar_id, gyration_colvar_id, mindist_colvar_id, &
        plane_distance_colvar_id, plane_plane_angle_colvar_id, population_colvar_id, &
        qparm_colvar_id, rmsd_colvar_id, torsion_colvar_id
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_compatible,&
                                              cp_unit_create,&
                                              cp_unit_from_cp2k,&
                                              cp_unit_release,&
                                              cp_unit_to_cp2k,&
                                              cp_unit_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_remove_values,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE input_val_types,                 ONLY: val_create,&
                                              val_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathlib,                         ONLY: invmat_symm
   USE md_environment_types,            ONLY: get_md_env,&
                                              md_environment_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                              fixd_constraint_type,&
                                              molecule_kind_type
   USE molecule_types,                  ONLY: global_constraint_type,&
                                              local_colvar_constraint_type,&
                                              molecule_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: massunit
   USE simpar_types,                    ONLY: simpar_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'blue_moon'
   PUBLIC :: blue_moon_initialize, blue_moon_update, blue_moon_write, blue_moon_save_restart, blue_moon_accumulate
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param md_env ...
!> \param constraint_section ...
!> \param start_step ...
!> \param restored ...
!> \param active ...
! **************************************************************************************************
   SUBROUTINE blue_moon_initialize(md_env, constraint_section, start_step, restored, active)
      TYPE(md_environment_type), POINTER :: md_env
      TYPE(section_vals_type), POINTER :: constraint_section
      INTEGER, INTENT(IN) :: start_step
      LOGICAL, INTENT(OUT), OPTIONAL :: restored
      LOGICAL, INTENT(OUT), OPTIONAL :: active
      TYPE(blue_moon_state_type), POINTER :: state
      TYPE(blue_moon_entry_type) :: tmp
      TYPE(force_env_type), POINTER :: force_env
      TYPE(cp_subsys_type), POINTER :: subsys
      TYPE(molecule_kind_type), POINTER :: kinds(:)
      TYPE(molecule_type), POINTER :: molecules(:)
      TYPE(global_constraint_type), POINTER :: gci
      TYPE(distribution_1d_type), POINTER :: local_molecules
      TYPE(mp_para_env_type), POINTER :: para_env
      TYPE(section_vals_type), POINTER :: restart, entries, ext_restart, print_key
      INTEGER :: n, i, j, k, m, nr, version, restart_step, cv_type
      INTEGER, POINTER :: ids(:), atoms(:)
      INTEGER, ALLOCATABLE :: owners(:)
      REAL(KIND=dp), POINTER :: values(:)
      CHARACTER(LEN=default_string_length) :: unit
      LOGICAL :: explicit, owned, ext_explicit, restart_counters, enabled, print_explicit

      CALL get_md_env(md_env, blue_moon=state, force_env=force_env, para_env=para_env)
      IF (PRESENT(restored)) restored = .FALSE.
      ! Explicit OFF print keys still request continuous work/restart bookkeeping.
      ! Ordinary constraint runs must not evaluate the additional joint metric.
      restart => section_vals_get_subs_vals(constraint_section, "BLUE_MOON_RESTART")
      CALL section_vals_get(restart, explicit=enabled)
      print_key => section_vals_get_subs_vals(constraint_section, "BLUE_MOON")
      CALL section_vals_get(print_key, explicit=print_explicit)
      enabled = enabled .OR. print_explicit
      print_key => section_vals_get_subs_vals(constraint_section, "COLVAR_SLOW_GROWTH")
      CALL section_vals_get(print_key, explicit=print_explicit)
      enabled = enabled .OR. print_explicit
      IF (PRESENT(active)) active = enabled
      CALL force_env_get(force_env, subsys=subsys)
      CALL cp_subsys_get(subsys, molecule_kind_set=kinds, molecule_set=molecules, gci=gci, &
                         local_molecules=local_molecules)
      gci%joint_colvar_checks = enabled
      IF (.NOT. enabled) RETURN
      n = 0
      DO m = 1, SIZE(molecules)
         IF (.NOT. ASSOCIATED(molecules(m)%lci)) CYCLE
         IF (.NOT. ASSOCIATED(molecules(m)%lci%lcolv)) CYCLE
         DO i = 1, SIZE(molecules(m)%lci%lcolv)
            IF (.NOT. molecules(m)%molecule_kind%colv_list(i)%restraint%active) n = n + 1
         END DO
      END DO
      IF (ASSOCIATED(gci%lcolv)) THEN
         DO i = 1, SIZE(gci%lcolv)
            IF (.NOT. gci%colv_list(i)%restraint%active) n = n + 1
         END DO
      END IF
      ALLOCATE (state%entries(n), owners(n))
      n = 0
      DO k = 1, SIZE(kinds)
         DO j = 1, kinds(k)%nmolecule
            m = kinds(k)%molecule_list(j)
            IF (.NOT. ASSOCIATED(molecules(m)%lci)) CYCLE
            IF (.NOT. ASSOCIATED(molecules(m)%lci%lcolv)) CYCLE
            owned = .FALSE.
            IF (local_molecules%n_el(k) > 0) owned = ANY(local_molecules%list(k)%array == m)
            CALL add_entries(molecules(m)%lci%lcolv, kinds(k)%colv_list, kinds(k)%fixd_list, m, owned)
         END DO
      END DO
      IF (ASSOCIATED(gci%lcolv)) THEN
         CALL add_entries(gci%lcolv, gci%colv_list, gci%fixd_list, 0, para_env%is_source())
      END IF
      ! Input order followed by global molecule ID is independent of MPI ownership.
      DO i = 2, n
         tmp = state%entries(i)
         j = i - 1
         DO WHILE (j >= 1)
            IF (.NOT. precedes(tmp, state%entries(j))) EXIT
            state%entries(j + 1) = state%entries(j)
            j = j - 1
         END DO
         state%entries(j + 1) = tmp
      END DO
      DO i = 1, n
         owners(i) = MERGE(1, 0, state%entries(i)%owned)
      END DO
      IF (n > 0) CALL para_env%sum(owners)
      IF (ANY(owners /= 1)) CPABORT("Blue Moon constraint must have exactly one MPI owner")
      state%step = start_step
      restart => section_vals_get_subs_vals(constraint_section, "BLUE_MOON_RESTART")
      CALL section_vals_get(restart, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(restart, "VERSION", i_val=version)
         CALL section_vals_val_get(restart, "STEP", i_val=restart_step)
         IF (version /= 1) CPABORT("Incompatible Blue Moon restart version")
         IF (restart_step /= start_step) THEN
            ext_restart => section_vals_get_subs_vals(force_env%root_section, "EXT_RESTART")
            CALL section_vals_get(ext_restart, explicit=ext_explicit)
            CALL section_vals_val_get(ext_restart, "RESTART_COUNTERS", l_val=restart_counters)
            IF (.NOT. ext_explicit .OR. restart_counters) THEN
               CPABORT("Incompatible Blue Moon restart MD step")
            END IF
         END IF
         entries => section_vals_get_subs_vals(restart, "CONSTRAINT")
         CALL section_vals_get(entries, n_repetition=nr)
         IF (nr /= n) CPABORT("Blue Moon restart constraint count changed")
         DO i = 1, n
            CALL section_vals_val_get(entries, "ID", i_vals=ids, i_rep_section=i)
            IF (SIZE(ids) /= 3) CPABORT("Invalid Blue Moon restart ID")
            IF (ANY(ids /= [state%entries(i)%input_index, state%entries(i)%molecule_index, &
                            state%entries(i)%local_index])) THEN
               CPABORT("Blue Moon restart constraint ID changed")
            END IF
            CALL section_vals_val_get(entries, "ATOMS", i_vals=atoms, i_rep_section=i)
            IF (SIZE(atoms) /= SIZE(state%entries(i)%local%colvar%i_atom)) THEN
               CPABORT("Blue Moon restart atom count changed")
            END IF
            IF (ANY(atoms /= state%entries(i)%local%colvar%i_atom)) THEN
               CPABORT("Blue Moon restart atoms changed")
            END IF
            CALL section_vals_val_get(entries, "CV_TYPE", i_val=cv_type, i_rep_section=i)
            IF (cv_type /= state%entries(i)%definition%type_id) THEN
               CPABORT("Blue Moon restart CV type changed")
            END IF
            CALL section_vals_val_get(entries, "CV_UNIT", c_val=unit, i_rep_section=i)
            CALL validate_cv_unit(state%entries(i), unit)
            state%entries(i)%cv_unit = unit
            state%entries(i)%definition%cv_unit = unit
            CALL section_vals_val_get(entries, "STATE", r_vals=values, i_rep_section=i)
            IF (SIZE(values) /= 9) CPABORT("Invalid Blue Moon restart state")
            state%entries(i)%previous_target = values(1)
            state%entries(i)%work = values(2)
            state%entries(i)%dwork = values(3)
            state%entries(i)%value = values(4)
            state%entries(i)%target = values(5)
            state%entries(i)%lambda_shake = values(6)
            state%entries(i)%lambda_rattle = values(7)
            state%entries(i)%gkt = values(8)
            state%entries(i)%target_velocity = values(9)
            IF (ABS(values(1) - state%entries(i)%definition%expected_value) > &
                1.0E-10_dp*MAX(1.0_dp, ABS(values(1)))) CPABORT("Blue Moon restart target changed")
            ! Official fixed-CV restarts initialize the integrator predictor to zero.
            ! Keep that behavior; moving constraints need the continuous pulling state.
            IF (state%entries(i)%definition%expected_value_growth_speed /= 0.0_dp) THEN
               state%entries(i)%local%lambda_shake = values(6)
               state%entries(i)%local%lambda_rattle = values(7)
               state%entries(i)%local%lambda = values(7)
            END IF
            state%entries(i)%definition%target_velocity = values(9)
         END DO
         CALL section_vals_val_get(restart, "LOG_WEIGHT", r_val=state%log_weight)
         state%restored = .TRUE.
      ELSE IF (start_step > 0 .AND. n > 0) THEN
         CPWARN("Old restart: Blue Moon work history unavailable; W starts a new segment")
      END IF
      state%initialized = .TRUE.
      IF (PRESENT(restored)) restored = state%restored
   CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param local ...
!> \param definitions ...
!> \param fixed ...
!> \param molecule_id ...
!> \param mine ...
! **************************************************************************************************
      SUBROUTINE add_entries(local, definitions, fixed, molecule_id, mine)
      TYPE(local_colvar_constraint_type), POINTER        :: local(:)
      TYPE(colvar_constraint_type), POINTER              :: definitions(:)
      TYPE(fixd_constraint_type), POINTER                :: fixed(:)
      INTEGER, INTENT(IN)                                :: molecule_id
      LOGICAL, INTENT(IN)                                :: mine

      INTEGER                                            :: a

         DO a = 1, SIZE(local)
            IF (definitions(a)%restraint%active) CYCLE
            IF (.NOT. ASSOCIATED(local(a)%colvar)) CPABORT("Missing active Blue Moon COLVAR")
            n = n + 1
            state%entries(n)%local => local(a)
            state%entries(n)%definition => definitions(a)
            state%entries(n)%fixed => fixed
            state%entries(n)%owned = mine
            state%entries(n)%input_index = definitions(a)%inp_seq_num
            state%entries(n)%molecule_index = molecule_id
            state%entries(n)%local_index = a
            state%entries(n)%cv_type = definitions(a)%type_id
            ALLOCATE (state%entries(n)%atoms(SIZE(local(a)%colvar%i_atom)))
            state%entries(n)%atoms(:) = local(a)%colvar%i_atom
            state%entries(n)%cv_unit = definitions(a)%cv_unit
            state%entries(n)%previous_target = definitions(a)%expected_value
            state%entries(n)%target = definitions(a)%expected_value
            state%entries(n)%value = local(a)%colvar%ss
         END DO
      END SUBROUTINE add_entries
   END SUBROUTINE blue_moon_initialize

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   FUNCTION precedes(a, b) RESULT(before)
      TYPE(blue_moon_entry_type), INTENT(IN) :: a, b
      INTEGER :: ai, bi
      LOGICAL :: before
      ai = MERGE(a%input_index, HUGE(0), a%input_index > 0)
      bi = MERGE(b%input_index, HUGE(0), b%input_index > 0)
      before = ai < bi
      IF (ai == bi) before = a%molecule_index < b%molecule_index
      IF (ai == bi .AND. a%molecule_index == b%molecule_index) before = a%local_index < b%local_index
   END FUNCTION precedes

! **************************************************************************************************
!> \brief ...
!> \param md_env ...
!> \param state ...
! **************************************************************************************************
   SUBROUTINE rebind_entries(md_env, state)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(blue_moon_state_type), INTENT(INOUT)          :: state

      INTEGER                                            :: i, j, k, m, n
      INTEGER, ALLOCATABLE                               :: owners(:)
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_type), POINTER                       :: molecules(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      ! Adaptive QM/MM can replace the subsystem, invalidating all borrowed pointers.
      CALL get_md_env(md_env, force_env=force_env, para_env=para_env)
      CALL force_env_get(force_env, subsys=subsys)
      CALL cp_subsys_get(subsys, molecule_set=molecules, gci=gci, local_molecules=local_molecules)
      n = 0
      DO m = 1, SIZE(molecules)
         IF (.NOT. ASSOCIATED(molecules(m)%lci)) CYCLE
         IF (.NOT. ASSOCIATED(molecules(m)%lci%lcolv)) CYCLE
         n = n + COUNT(.NOT. molecules(m)%molecule_kind%colv_list%restraint%active)
      END DO
      IF (ASSOCIATED(gci%lcolv)) n = n + COUNT(.NOT. gci%colv_list%restraint%active)
      IF (n /= SIZE(state%entries)) CPABORT("Blue Moon constraint topology changed during MD")
      ALLOCATE (owners(n))
      DO i = 1, n
         m = state%entries(i)%molecule_index
         j = state%entries(i)%local_index
         IF (m == 0) THEN
            IF (.NOT. ASSOCIATED(gci%lcolv)) CPABORT("Blue Moon global constraints disappeared")
            IF (j < 1 .OR. j > SIZE(gci%lcolv)) CPABORT("Blue Moon constraint index changed")
            state%entries(i)%local => gci%lcolv(j)
            state%entries(i)%definition => gci%colv_list(j)
            state%entries(i)%fixed => gci%fixd_list
            state%entries(i)%owned = para_env%is_source()
         ELSE
            IF (m < 1 .OR. m > SIZE(molecules)) CPABORT("Blue Moon molecule index changed")
            IF (.NOT. ASSOCIATED(molecules(m)%lci)) CPABORT("Blue Moon local constraints disappeared")
            IF (.NOT. ASSOCIATED(molecules(m)%lci%lcolv)) THEN
               CPABORT("Blue Moon local constraints disappeared")
            END IF
            IF (j < 1 .OR. j > SIZE(molecules(m)%lci%lcolv)) THEN
               CPABORT("Blue Moon constraint index changed")
            END IF
            state%entries(i)%local => molecules(m)%lci%lcolv(j)
            state%entries(i)%definition => molecules(m)%molecule_kind%colv_list(j)
            state%entries(i)%fixed => molecules(m)%molecule_kind%fixd_list
            k = molecules(m)%molecule_kind%kind_number
            state%entries(i)%owned = .FALSE.
            IF (local_molecules%n_el(k) > 0) THEN
               state%entries(i)%owned = ANY(local_molecules%list(k)%array == m)
            END IF
         END IF
         IF (state%entries(i)%definition%restraint%active .OR. &
             state%entries(i)%definition%inp_seq_num /= state%entries(i)%input_index .OR. &
             state%entries(i)%definition%type_id /= state%entries(i)%cv_type) THEN
            CPABORT("Blue Moon constraint identity changed during MD")
         END IF
         IF (SIZE(state%entries(i)%local%colvar%i_atom) /= SIZE(state%entries(i)%atoms)) THEN
            CPABORT("Blue Moon constraint atoms changed during MD")
         END IF
         IF (ANY(state%entries(i)%local%colvar%i_atom /= state%entries(i)%atoms)) THEN
            CPABORT("Blue Moon constraint atoms changed during MD")
         END IF
         state%entries(i)%definition%cv_unit = state%entries(i)%cv_unit
         owners(i) = MERGE(1, 0, state%entries(i)%owned)
      END DO
      IF (n > 0) CALL para_env%sum(owners)
      IF (ANY(owners /= 1)) CPABORT("Blue Moon constraint must have exactly one MPI owner")
   END SUBROUTINE rebind_entries

! **************************************************************************************************
!> \brief ...
!> \param md_env ...
! **************************************************************************************************
   SUBROUTINE blue_moon_accumulate(md_env)
      TYPE(md_environment_type), POINTER                 :: md_env

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: increment
      REAL(KIND=dp), ALLOCATABLE                         :: lambdas(:)
      TYPE(blue_moon_state_type), POINTER                :: state
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL get_md_env(md_env, blue_moon=state, para_env=para_env)
      IF (.NOT. state%initialized) RETURN
      n = SIZE(state%entries)
      IF (n == 0) RETURN
      CALL rebind_entries(md_env, state)
      ALLOCATE (lambdas(n))
      lambdas = 0.0_dp
      DO i = 1, n
         IF (state%entries(i)%owned) lambdas(i) = state%entries(i)%local%lambda_shake
      END DO
      CALL para_env%sum(lambdas)
      DO i = 1, n
         IF (.NOT. state%work_pending) state%entries(i)%dwork = 0.0_dp
         increment = -lambdas(i)*(state%entries(i)%definition%expected_value - state%entries(i)%previous_target)
         state%entries(i)%dwork = state%entries(i)%dwork + increment
         state%entries(i)%work = state%entries(i)%work + increment
         state%entries(i)%previous_target = state%entries(i)%definition%expected_value
      END DO
      state%work_pending = .TRUE.
   END SUBROUTINE blue_moon_accumulate

! **************************************************************************************************
!> \brief ...
!> \param md_env ...
!> \param step ...
! **************************************************************************************************
   SUBROUTINE blue_moon_update(md_env, step)
      TYPE(md_environment_type), POINTER                 :: md_env
      INTEGER, INTENT(IN)                                :: step

      INTEGER                                            :: i, n
      REAL(KIND=dp), ALLOCATABLE                         :: correction(:), multipliers(:, :)
      TYPE(blue_moon_state_type), POINTER                :: state
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(simpar_type), POINTER                         :: simpar

      CALL get_md_env(md_env, blue_moon=state, force_env=force_env, para_env=para_env, simpar=simpar)
      IF (.NOT. state%initialized) RETURN
      n = SIZE(state%entries)
      IF (n == 0) THEN
         state%step = step
         RETURN
      END IF
      IF (step == state%step) RETURN
      IF (step /= state%step + 1) CPABORT("Blue Moon state must be updated once per MD step")
      CALL rebind_entries(md_env, state)
      IF (.NOT. state%work_pending) CALL blue_moon_accumulate(md_env)
      ALLOCATE (multipliers(2, n), correction(n))
      multipliers = 0.0_dp
      DO i = 1, n
         IF (.NOT. state%entries(i)%owned) CYCLE
         multipliers(1, i) = state%entries(i)%local%lambda_shake
         multipliers(2, i) = state%entries(i)%local%lambda_rattle
      END DO
      IF (n > 0) CALL para_env%sum(multipliers)
      CALL force_env_get(force_env, subsys=subsys)
      CALL cp_subsys_get(subsys, particle_set=particles, cell=cell)
      IF (n > 0) CALL joint_metric(state, particles, cell, simpar%temp_ext, state%log_weight, correction)
      DO i = 1, n
         state%entries(i)%lambda_shake = multipliers(1, i)
         state%entries(i)%lambda_rattle = multipliers(2, i)
         state%entries(i)%value = state%entries(i)%local%colvar%ss
         state%entries(i)%target = state%entries(i)%definition%expected_value
         state%entries(i)%target_velocity = state%entries(i)%definition%target_velocity
         state%entries(i)%gkt = correction(i)
      END DO
      state%step = step
      state%work_pending = .FALSE.
   END SUBROUTINE blue_moon_update

! **************************************************************************************************
!> \brief ...
!> \param state ...
!> \param particles ...
!> \param cell ...
!> \param pos ...
!> \param bmat ...
! **************************************************************************************************
   SUBROUTINE evaluate(state, particles, cell, pos, bmat)
      TYPE(blue_moon_state_type), INTENT(IN)             :: state
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: pos(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: bmat(:, :, :)

      INTEGER                                            :: a, i, j
      REAL(KIND=dp)                                      :: angstrom_scale

      angstrom_scale = cp_unit_from_cp2k(1.0_dp, "angstrom")
      bmat = 0.0_dp
      DO i = 1, SIZE(state%entries)
         CALL colvar_eval_mol_f(state%entries(i)%local%colvar, cell, particles=particles, pos=pos, &
                                fixd_list=state%entries(i)%fixed)
         DO j = 1, SIZE(state%entries(i)%local%colvar%i_atom)
            a = state%entries(i)%local%colvar%i_atom(j)
            IF (a < 1 .OR. a > SIZE(particles)) CPABORT("Invalid Blue Moon COLVAR atom index")
            bmat(i, :, a) = bmat(i, :, a) + state%entries(i)%local%colvar%dsdr(:, j)/angstrom_scale
         END DO
      END DO
   END SUBROUTINE evaluate

! **************************************************************************************************
!> \brief ...
!> \param bmat ...
!> \param masses ...
!> \param inverse ...
!> \param logdet ...
! **************************************************************************************************
   SUBROUTINE metric(bmat, masses, inverse, logdet)
      REAL(KIND=dp), INTENT(IN)                          :: bmat(:, :, :), masses(:)
      REAL(KIND=dp), INTENT(OUT)                         :: inverse(:, :), logdet

      INTEGER                                            :: a, i, info, j, n
      INTEGER, ALLOCATABLE                               :: iwork(:)
      REAL(KIND=dp)                                      :: norm, rcond
      REAL(KIND=dp), ALLOCATABLE                         :: scales(:), work(:), z(:, :)

      n = SIZE(bmat, 1)
      ALLOCATE (z(n, n), scales(n), work(3*n), iwork(n))
      z = 0.0_dp
      DO a = 1, SIZE(masses)
         IF (masses(a) <= 0.0_dp) CPABORT("Blue Moon requires positive atomic masses")
         DO i = 1, n
            DO j = 1, n
               z(i, j) = z(i, j) + DOT_PRODUCT(bmat(i, :, a), bmat(j, :, a))/masses(a)
            END DO
         END DO
      END DO
      DO i = 1, n
         IF (z(i, i) <= 0.0_dp) CPABORT("Blue Moon CV has a zero mass-metric gradient")
         scales(i) = SQRT(z(i, i))
      END DO
      DO i = 1, n
         DO j = 1, n
            inverse(i, j) = z(i, j)/scales(i)/scales(j)
         END DO
      END DO
      norm = MAXVAL(SUM(ABS(inverse), DIM=1))
      CALL dpotrf('U', n, inverse, n, info)
      IF (info /= 0) CPABORT("Dependent Blue Moon constraints: mass metric is not positive definite")
      CALL dpocon('U', n, inverse, n, norm, rcond, work, iwork, info)
      IF (info /= 0 .OR. rcond < 100.0_dp*EPSILON(1.0_dp)) THEN
         CPABORT("Dependent or ill-conditioned Blue Moon constraints")
      END IF
      logdet = 2.0_dp*SUM(LOG(scales))
      DO i = 1, n
         logdet = logdet + 2.0_dp*LOG(inverse(i, i))
      END DO
      CALL invmat_symm(inverse, potrf=.FALSE., uplo='U')
      DO i = 1, n
         DO j = 1, n
            inverse(i, j) = inverse(i, j)/scales(i)/scales(j)
         END DO
      END DO
   END SUBROUTINE metric

! **************************************************************************************************
!> \brief ...
!> \param state ...
!> \param particles ...
!> \param cell ...
!> \param kbt ...
!> \param log_weight ...
!> \param gkt ...
! **************************************************************************************************
   SUBROUTINE joint_metric(state, particles, cell, kbt, log_weight, gkt)
      TYPE(blue_moon_state_type), INTENT(IN)             :: state
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: kbt
      REAL(KIND=dp), INTENT(OUT)                         :: log_weight, gkt(:)

      INTEGER                                            :: a, d, n, na
      REAL(KIND=dp)                                      :: derivative, h, logdet, minus, orig, plus
      REAL(KIND=dp), ALLOCATABLE                         :: b(:, :, :), contraction(:), &
                                                            inverse(:, :), masses(:), pos(:, :), &
                                                            trial(:, :, :), trial_inverse(:, :)

      n = SIZE(state%entries)
      na = SIZE(particles)
      ALLOCATE (pos(3, na), b(n, 3, na), trial(n, 3, na), masses(na), &
                inverse(n, n), trial_inverse(n, n), contraction(n))
      DO a = 1, na
         pos(:, a) = particles(a)%r
         CALL get_atomic_kind(particles(a)%atomic_kind, mass=masses(a))
      END DO
      masses = masses/massunit
      CALL evaluate(state, particles, cell, pos, b)
      CALL metric(b, masses, inverse, logdet)
      log_weight = -0.5_dp*logdet
      contraction = 0.0_dp
      h = cp_unit_to_cp2k(1.0E-4_dp, "angstrom")
      DO a = 1, na
         IF (ALL(b(:, :, a) == 0.0_dp)) CYCLE
         DO d = 1, 3
            orig = pos(d, a)
            pos(d, a) = orig + h
            CALL evaluate(state, particles, cell, pos, trial)
            CALL metric(trial, masses, trial_inverse, plus)
            pos(d, a) = orig - h
            CALL evaluate(state, particles, cell, pos, trial)
            CALL metric(trial, masses, trial_inverse, minus)
            pos(d, a) = orig
            derivative = (plus - minus)/(2.0E-4_dp)
            contraction(:) = contraction + b(:, d, a)*derivative/masses(a)
         END DO
      END DO
      CALL evaluate(state, particles, cell, pos, trial)
      ! Trial evaluations modify both particle positions and the CV derivative cache.
      DO a = 1, na
         particles(a)%r = pos(:, a)
      END DO
      gkt = 0.5_dp*kbt*MATMUL(inverse, contraction)
   END SUBROUTINE joint_metric

! **************************************************************************************************
!> \brief ...
!> \param ENTRY ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_cv_label(ENTRY) RESULT(label)
      TYPE(blue_moon_entry_type), INTENT(IN)             :: entry
      CHARACTER(LEN=default_string_length)               :: label

      SELECT CASE (entry%definition%type_id)
      CASE (dist_colvar_id, plane_distance_colvar_id, gyration_colvar_id, rmsd_colvar_id, mindist_colvar_id)
         label = "bohr"
      CASE (angle_colvar_id, torsion_colvar_id, plane_plane_angle_colvar_id)
         label = "rad"
      CASE (coord_colvar_id, population_colvar_id, qparm_colvar_id)
         label = "1"
      CASE DEFAULT
         label = "internal_cp2k"
      END SELECT
   END FUNCTION native_cv_label

! **************************************************************************************************
!> \brief ...
!> \param ENTRY ...
!> \param unit ...
! **************************************************************************************************
   SUBROUTINE validate_cv_unit(ENTRY, unit)
      TYPE(blue_moon_entry_type), INTENT(IN)             :: entry
      CHARACTER(LEN=*), INTENT(IN)                       :: unit

      CHARACTER(LEN=default_string_length)               :: reference
      TYPE(cp_unit_type)                                 :: candidate, ref

      IF (unit == "internal_cp2k") RETURN
      reference = entry%cv_unit
      IF (reference == "internal_cp2k") reference = native_cv_label(ENTRY)
      IF (reference == "internal_cp2k") THEN
         CPABORT("CV_UNIT override requires a dimensioned input TARGET for this COLVAR")
      END IF
      IF (reference == "1") reference = ""
      CALL cp_unit_create(ref, TRIM(reference))
      IF (unit == "1") THEN
         CALL cp_unit_create(candidate, "")
      ELSE
         CALL cp_unit_create(candidate, TRIM(unit))
      END IF
      IF (.NOT. cp_unit_compatible(ref, candidate) .OR. ANY(ref%power /= candidate%power)) THEN
         CPABORT("CV_UNIT is incompatible with a Blue Moon COLVAR")
      END IF
      CALL cp_unit_release(ref)
      CALL cp_unit_release(candidate)
   END SUBROUTINE validate_cv_unit

! **************************************************************************************************
!> \brief ...
!> \param reference ...
!> \param unit ...
! **************************************************************************************************
   SUBROUTINE validate_display_unit(reference, unit)
      CHARACTER(LEN=*), INTENT(IN)                       :: reference, unit

      TYPE(cp_unit_type)                                 :: candidate, ref

      IF (unit == "internal_cp2k") RETURN
      CALL cp_unit_create(ref, reference)
      CALL cp_unit_create(candidate, unit)
      IF (.NOT. cp_unit_compatible(ref, candidate) .OR. ANY(ref%power /= candidate%power)) THEN
         CPABORT("Invalid Blue Moon ENERGY_UNIT or TIME_UNIT dimension")
      END IF
      CALL cp_unit_release(ref)
      CALL cp_unit_release(candidate)
   END SUBROUTINE validate_display_unit

! **************************************************************************************************
!> \brief ...
!> \param md_env ...
!> \param constraint_section ...
!> \param print_key_name ...
!> \param unit_nr ...
!> \param is_new_file ...
!> \param time ...
!> \param step ...
! **************************************************************************************************
   SUBROUTINE blue_moon_write(md_env, constraint_section, print_key_name, unit_nr, is_new_file, time, step)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(section_vals_type), POINTER                   :: constraint_section
      CHARACTER(LEN=*), INTENT(IN)                       :: print_key_name
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: is_new_file
      REAL(KIND=dp), INTENT(IN)                          :: time
      INTEGER, INTENT(IN)                                :: step

      CHARACTER(LEN=16)                                  :: count
      CHARACTER(LEN=:), ALLOCATABLE                      :: metric_unit, product
      CHARACTER(LEN=default_string_length)               :: energy_label, energy_unit, requested, &
                                                            time_label, time_unit
      CHARACTER(LEN=default_string_length), ALLOCATABLE  :: forces(:), labels(:), units(:)
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: energy_scale, force, gkt, log_weight, &
                                                            TARGET, value, weight
      REAL(KIND=dp), ALLOCATABLE                         :: cv_scale(:)
      TYPE(blue_moon_state_type), POINTER                :: state
      TYPE(section_vals_type), POINTER                   :: print_key

      IF (unit_nr <= 0) RETURN
      CALL get_md_env(md_env, blue_moon=state)
      IF (.NOT. state%initialized) RETURN
      CPASSERT(state%step == step)
      n = SIZE(state%entries)
      IF (n == 0) RETURN
      print_key => section_vals_get_subs_vals(constraint_section, print_key_name)
      CALL section_vals_val_get(print_key, "CV_UNIT", c_val=requested)
      CALL section_vals_val_get(print_key, "ENERGY_UNIT", c_val=energy_unit)
      CALL section_vals_val_get(print_key, "TIME_UNIT", c_val=time_unit)
      time_label = time_unit
      IF (time_label == "internal_cp2k") time_label = "au_t"
      CALL validate_display_unit("hartree", TRIM(energy_unit))
      CALL validate_display_unit("fs", TRIM(time_unit))
      ALLOCATE (units(n), labels(n), forces(n), cv_scale(n))
      energy_scale = cp_unit_from_cp2k(1.0_dp, TRIM(energy_unit))
      energy_label = energy_unit
      IF (energy_label == "internal_cp2k") energy_label = "hartree"
      log_weight = state%log_weight
      product = ""
      DO i = 1, n
         units(i) = state%entries(i)%cv_unit
         IF (requested /= "input") THEN
            CALL validate_cv_unit(state%entries(i), TRIM(requested))
            units(i) = requested
         END IF
         labels(i) = units(i)
         IF (labels(i) == "internal_cp2k") labels(i) = native_cv_label(state%entries(i))
         forces(i) = TRIM(energy_label)//"*"//TRIM(labels(i))//"^-1"
         cv_scale(i) = 1.0_dp
         IF (units(i) /= "1") cv_scale(i) = cp_unit_from_cp2k(1.0_dp, TRIM(units(i)))
         log_weight = log_weight - LOG(cv_scale(i))
         IF (i > 1) product = product//"*"
         product = product//TRIM(labels(i))
      END DO
      IF (log_weight > LOG(HUGE(1.0_dp)) .OR. log_weight < LOG(TINY(1.0_dp))) THEN
         CPABORT("Blue Moon output metric weight is outside representable range")
      END IF
      weight = EXP(log_weight)
      WRITE (count, '(I0)') n
      metric_unit = "angstrom^"//TRIM(count)//"*amu^("//TRIM(count)//"/2)/("//product//")"
      IF (n == 1) metric_unit = "angstrom*amu^0.5*"//TRIM(labels(1))//"^-1"
      IF (is_new_file) THEN
         WRITE (unit_nr, '(A)') "# CP2K Blue Moon/TI COLVAR output"
         WRITE (unit_nr, '(A)') "# TI estimator = <blue_moon_force>/<inv_sqrt_z>; W is mechanical pulling work"
         WRITE (unit_nr, '(A)') "# z convention = all active CV gradients per R[angstrom], mass[amu], joint determinant"
         WRITE (unit_nr, '(A,A)') "# time_unit = ", TRIM(time_label)
         WRITE (unit_nr, '(A,A)') "# energy_unit = ", TRIM(energy_label)
         WRITE (unit_nr, '(A,A)') "# metric_unit = ", metric_unit
         IF (n == 1) THEN
            WRITE (unit_nr, '(A,A)') "# cv_unit = ", TRIM(units(1))
            WRITE (unit_nr, '(A,A)') "# native_cv_unit = ", TRIM(native_cv_label(state%entries(1)))
            WRITE (unit_nr, '(A,A)') "# force_unit = ", TRIM(forces(1))
            WRITE (unit_nr, '(A,A)') "# blue_moon_force_unit = ", TRIM(forces(1))//"*("//metric_unit//")"
         ELSE
            WRITE (unit_nr, '(A)') "# cv_unit = per_constraint"
            WRITE (unit_nr, '(A)') "# force_unit = per_constraint"
         END IF
         DO i = 1, n
            WRITE (unit_nr, '(A,I0,A,3(I0,1X))') "# constraint ", i, " id = ", &
               state%entries(i)%input_index, state%entries(i)%molecule_index, state%entries(i)%local_index
            WRITE (unit_nr, '(A,I0,A,A)') "# constraint ", i, " cv_unit = ", TRIM(units(i))
            WRITE (unit_nr, '(A,I0,A,A)') "# constraint ", i, " native_cv_unit = ", TRIM(native_cv_label(state%entries(i)))
            WRITE (unit_nr, '(A,I0,A,A)') "# constraint ", i, " force_unit = ", TRIM(forces(i))
            WRITE (unit_nr, '(A,I0,A,A)') "# constraint ", i, " blue_moon_force_unit = ", &
               TRIM(forces(i))//"*("//metric_unit//")"
         END DO
         WRITE (unit_nr, '(A)') "# columns: step time constraint value target lambda inv_sqrt_z gkt blue_moon_force dxi dW W"
      END IF
      DO i = 1, n
         value = state%entries(i)%value*cv_scale(i)
         TARGET = state%entries(i)%target*cv_scale(i)
         force = -state%entries(i)%lambda_shake*energy_scale/cv_scale(i)
         gkt = state%entries(i)%gkt*energy_scale/cv_scale(i)
         WRITE (unit_nr, '(I10,1X,ES20.10,1X,I8,1X,9(ES24.14,1X))') &
            step, cp_unit_from_cp2k(time, TRIM(time_unit)), i, value, TARGET, force, weight, gkt, &
            weight*(force + gkt), value - TARGET, state%entries(i)%dwork*energy_scale, state%entries(i)%work*energy_scale
      END DO
   END SUBROUTINE blue_moon_write

! **************************************************************************************************
!> \brief ...
!> \param md_env ...
!> \param motion_section ...
! **************************************************************************************************
   SUBROUTINE blue_moon_save_restart(md_env, motion_section)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(section_vals_type), POINTER                   :: motion_section

      INTEGER                                            :: i
      TYPE(blue_moon_state_type), POINTER                :: state
      TYPE(section_vals_type), POINTER                   :: entries, restart
      TYPE(val_type), POINTER                            :: val

      CALL get_md_env(md_env, blue_moon=state)
      IF (.NOT. state%initialized) RETURN
      IF (SIZE(state%entries) == 0) RETURN
      restart => section_vals_get_subs_vals(motion_section, "CONSTRAINT%BLUE_MOON_RESTART")
      CALL section_vals_remove_values(restart)
      CALL section_vals_val_set(restart, "VERSION", i_val=1)
      CALL section_vals_val_set(restart, "STEP", i_val=state%step)
      CALL section_vals_val_set(restart, "LOG_WEIGHT", r_val=state%log_weight)
      entries => section_vals_get_subs_vals(restart, "CONSTRAINT")
      NULLIFY (val)
      DO i = 1, SIZE(state%entries)
         CALL val_create(val, i_vals=[state%entries(i)%input_index, state%entries(i)%molecule_index, &
                                      state%entries(i)%local_index])
         CALL section_vals_val_set(entries, "ID", val=val, i_rep_section=i)
         NULLIFY (val)
         CALL val_create(val, i_vals=state%entries(i)%local%colvar%i_atom)
         CALL section_vals_val_set(entries, "ATOMS", val=val, i_rep_section=i)
         NULLIFY (val)
         CALL section_vals_val_set(entries, "CV_TYPE", i_val=state%entries(i)%definition%type_id, i_rep_section=i)
         CALL section_vals_val_set(entries, "CV_UNIT", c_val=TRIM(state%entries(i)%cv_unit), i_rep_section=i)
         CALL val_create(val, r_vals=[state%entries(i)%previous_target, state%entries(i)%work, &
                                      state%entries(i)%dwork, state%entries(i)%value, state%entries(i)%target, &
                                      state%entries(i)%lambda_shake, state%entries(i)%lambda_rattle, &
                                      state%entries(i)%gkt, state%entries(i)%target_velocity])
         CALL section_vals_val_set(entries, "STATE", val=val, i_rep_section=i)
         NULLIFY (val)
      END DO
   END SUBROUTINE blue_moon_save_restart
END MODULE blue_moon
