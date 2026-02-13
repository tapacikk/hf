CC      = cc
CFLAGS  = -O2 -fopenmp -I/opt/libcint/latest/include
LDFLAGS = -L/opt/libcint/latest/lib -Wl,-rpath=/opt/libcint/latest/lib
LDLIBS  = -lcint

TARGET  = i_wat_dz
SRC     = i_wat_dz.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET)
