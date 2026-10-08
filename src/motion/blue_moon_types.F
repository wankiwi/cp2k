!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE blue_moon_types
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                              fixd_constraint_type
   USE molecule_types,                  ONLY: local_colvar_constraint_type

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: blue_moon_entry_type, blue_moon_state_type

   TYPE blue_moon_entry_type
      INTEGER :: input_index = 0, molecule_index = 0, local_index = 0
      INTEGER :: cv_type = 0
      INTEGER, ALLOCATABLE :: atoms(:)
      LOGICAL :: owned = .FALSE.
      CHARACTER(LEN=default_string_length) :: cv_unit = "internal_cp2k"
      TYPE(local_colvar_constraint_type), POINTER :: local => NULL()
      TYPE(colvar_constraint_type), POINTER :: definition => NULL()
      TYPE(fixd_constraint_type), POINTER :: fixed(:) => NULL()
      REAL(KIND=dp) :: previous_target = 0.0_dp, work = 0.0_dp, dwork = 0.0_dp
      REAL(KIND=dp) :: value = 0.0_dp, TARGET = 0.0_dp, lambda_shake = 0.0_dp
      REAL(KIND=dp) :: lambda_rattle = 0.0_dp, gkt = 0.0_dp, target_velocity = 0.0_dp
   END TYPE blue_moon_entry_type

   TYPE blue_moon_state_type
      LOGICAL :: initialized = .FALSE., restored = .FALSE., work_pending = .FALSE.
      INTEGER :: step = -1
      REAL(KIND=dp) :: log_weight = 0.0_dp
      TYPE(blue_moon_entry_type), ALLOCATABLE :: entries(:)
   END TYPE blue_moon_state_type
END MODULE blue_moon_types
