FC      := gfortran
FFLAGS  := -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none
LDFLAGS :=

SRC_DIR := src
OBJ_DIR := build
TARGET  := m_poisson_mg.x

SRCS := \
	$(SRC_DIR)/precision.f90 \
	$(SRC_DIR)/m_boundary_condition.f90 \
	$(SRC_DIR)/m_poisson_bc_ops.f90 \
	$(SRC_DIR)/m_grid_obj.f90 \
	$(SRC_DIR)/m_poison_mg.f90 \
	$(SRC_DIR)/m_iorho_serial.f90 \
	$(SRC_DIR)/m_poisson_mg.F

OBJS := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(filter %.f90,$(SRCS))) \
        $(patsubst $(SRC_DIR)/%.F,$(OBJ_DIR)/%.o,$(filter %.F,$(SRCS)))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(OBJS) -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -J$(OBJ_DIR) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F | $(OBJ_DIR)
	$(FC) $(FFLAGS) -ffree-form -c $< -J$(OBJ_DIR) -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean
