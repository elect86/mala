// Generated by jextract

package lammps;

import java.lang.invoke.*;
import java.lang.foreign.*;
import java.nio.ByteOrder;
import java.nio.file.Path;
import java.util.*;
import java.util.function.*;
import java.util.stream.*;

import static java.lang.foreign.ValueLayout.*;
import static java.lang.foreign.MemoryLayout.PathElement.*;

public class library_h {

    library_h() {
        // Should not be called directly
    }

    static final Arena LIBRARY_ARENA = Arena.ofAuto();
    static final boolean TRACE_DOWNCALLS = Boolean.getBoolean("jextract.trace.downcalls");

    static void traceDowncall(String name, Object... args) {
         String traceArgs = Arrays.stream(args)
                       .map(Object::toString)
                       .collect(Collectors.joining(", "));
         System.out.printf("%s(%s)\n", name, traceArgs);
    }

    static MemorySegment findOrThrow(String symbol) {
        return SYMBOL_LOOKUP.find(symbol)
            .orElseThrow(() -> new UnsatisfiedLinkError("unresolved symbol: " + symbol));
    }

    static MethodHandle upcallHandle(Class<?> fi, String name, FunctionDescriptor fdesc) {
        try {
            return MethodHandles.lookup().findVirtual(fi, name, fdesc.toMethodType());
        } catch (ReflectiveOperationException ex) {
            throw new AssertionError(ex);
        }
    }

    static MemoryLayout align(MemoryLayout layout, long align) {
        return switch (layout) {
            case PaddingLayout p -> p;
            case ValueLayout v -> v.withByteAlignment(align);
            case GroupLayout g -> {
                MemoryLayout[] alignedMembers = g.memberLayouts().stream()
                        .map(m -> align(m, align)).toArray(MemoryLayout[]::new);
                yield g instanceof StructLayout ?
                        MemoryLayout.structLayout(alignedMembers) : MemoryLayout.unionLayout(alignedMembers);
            }
            case SequenceLayout s -> MemoryLayout.sequenceLayout(s.elementCount(), align(s.elementLayout(), align));
        };
    }

    static final SymbolLookup SYMBOL_LOOKUP = SymbolLookup.libraryLookup(Path.of("/home/elect/PycharmProjects/lammps/build/liblammps.so"), LIBRARY_ARENA)
            .or(SymbolLookup.loaderLookup())
            .or(Linker.nativeLinker().defaultLookup());

    public static final ValueLayout.OfBoolean C_BOOL = ValueLayout.JAVA_BOOLEAN;
    public static final ValueLayout.OfByte C_CHAR = ValueLayout.JAVA_BYTE;
    public static final ValueLayout.OfShort C_SHORT = ValueLayout.JAVA_SHORT;
    public static final ValueLayout.OfInt C_INT = ValueLayout.JAVA_INT;
    public static final ValueLayout.OfLong C_LONG_LONG = ValueLayout.JAVA_LONG;
    public static final ValueLayout.OfFloat C_FLOAT = ValueLayout.JAVA_FLOAT;
    public static final ValueLayout.OfDouble C_DOUBLE = ValueLayout.JAVA_DOUBLE;
    public static final AddressLayout C_POINTER = ValueLayout.ADDRESS
            .withTargetLayout(MemoryLayout.sequenceLayout(java.lang.Long.MAX_VALUE, JAVA_BYTE));
    public static final ValueLayout.OfLong C_LONG = ValueLayout.JAVA_LONG;

    private static class lammps_open_no_mpi {
        public static final FunctionDescriptor DESC = FunctionDescriptor.of(
            library_h.C_POINTER,
            library_h.C_INT,
            library_h.C_POINTER,
            library_h.C_POINTER
        );

        public static final MemorySegment ADDR = library_h.findOrThrow("lammps_open_no_mpi");

        public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
    }

    /**
     * Function descriptor for:
     * {@snippet lang=c :
     * void *lammps_open_no_mpi(int argc, char **argv, void **ptr)
     * }
     */
    public static FunctionDescriptor lammps_open_no_mpi$descriptor() {
        return lammps_open_no_mpi.DESC;
    }

    /**
     * Downcall method handle for:
     * {@snippet lang=c :
     * void *lammps_open_no_mpi(int argc, char **argv, void **ptr)
     * }
     */
    public static MethodHandle lammps_open_no_mpi$handle() {
        return lammps_open_no_mpi.HANDLE;
    }

    /**
     * Address for:
     * {@snippet lang=c :
     * void *lammps_open_no_mpi(int argc, char **argv, void **ptr)
     * }
     */
    public static MemorySegment lammps_open_no_mpi$address() {
        return lammps_open_no_mpi.ADDR;
    }

    /**
     * {@snippet lang=c :
     * void *lammps_open_no_mpi(int argc, char **argv, void **ptr)
     * }
     */
    public static MemorySegment lammps_open_no_mpi(int argc, MemorySegment argv, MemorySegment ptr) {
        var mh$ = lammps_open_no_mpi.HANDLE;
        try {
            if (TRACE_DOWNCALLS) {
                traceDowncall("lammps_open_no_mpi", argc, argv, ptr);
            }
            return (MemorySegment)mh$.invokeExact(argc, argv, ptr);
        } catch (Throwable ex$) {
           throw new AssertionError("should not reach here", ex$);
        }
    }

    private static class lammps_file {
        public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
            library_h.C_POINTER,
            library_h.C_POINTER
        );

        public static final MemorySegment ADDR = library_h.findOrThrow("lammps_file");

        public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
    }

    /**
     * Function descriptor for:
     * {@snippet lang=c :
     * void lammps_file(void *handle, const char *file)
     * }
     */
    public static FunctionDescriptor lammps_file$descriptor() {
        return lammps_file.DESC;
    }

    /**
     * Downcall method handle for:
     * {@snippet lang=c :
     * void lammps_file(void *handle, const char *file)
     * }
     */
    public static MethodHandle lammps_file$handle() {
        return lammps_file.HANDLE;
    }

    /**
     * Address for:
     * {@snippet lang=c :
     * void lammps_file(void *handle, const char *file)
     * }
     */
    public static MemorySegment lammps_file$address() {
        return lammps_file.ADDR;
    }

    /**
     * {@snippet lang=c :
     * void lammps_file(void *handle, const char *file)
     * }
     */
    public static void lammps_file(MemorySegment handle, MemorySegment file) {
        var mh$ = lammps_file.HANDLE;
        try {
            if (TRACE_DOWNCALLS) {
                traceDowncall("lammps_file", handle, file);
            }
            mh$.invokeExact(handle, file);
        } catch (Throwable ex$) {
           throw new AssertionError("should not reach here", ex$);
        }
    }

    private static class lammps_extract_compute {
        public static final FunctionDescriptor DESC = FunctionDescriptor.of(
            library_h.C_POINTER,
            library_h.C_POINTER,
            library_h.C_POINTER,
            library_h.C_INT,
            library_h.C_INT
        );

        public static final MemorySegment ADDR = library_h.findOrThrow("lammps_extract_compute");

        public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
    }

    /**
     * Function descriptor for:
     * {@snippet lang=c :
     * void *lammps_extract_compute(void *handle, const char *, int, int)
     * }
     */
    public static FunctionDescriptor lammps_extract_compute$descriptor() {
        return lammps_extract_compute.DESC;
    }

    /**
     * Downcall method handle for:
     * {@snippet lang=c :
     * void *lammps_extract_compute(void *handle, const char *, int, int)
     * }
     */
    public static MethodHandle lammps_extract_compute$handle() {
        return lammps_extract_compute.HANDLE;
    }

    /**
     * Address for:
     * {@snippet lang=c :
     * void *lammps_extract_compute(void *handle, const char *, int, int)
     * }
     */
    public static MemorySegment lammps_extract_compute$address() {
        return lammps_extract_compute.ADDR;
    }

    /**
     * {@snippet lang=c :
     * void *lammps_extract_compute(void *handle, const char *, int, int)
     * }
     */
    public static MemorySegment lammps_extract_compute(MemorySegment handle, MemorySegment x1, int x2, int x3) {
        var mh$ = lammps_extract_compute.HANDLE;
        try {
            if (TRACE_DOWNCALLS) {
                traceDowncall("lammps_extract_compute", handle, x1, x2, x3);
            }
            return (MemorySegment)mh$.invokeExact(handle, x1, x2, x3);
        } catch (Throwable ex$) {
           throw new AssertionError("should not reach here", ex$);
        }
    }

    public enum Style { global, atom, local }
    public enum Type { scalar, vector, array, sizeVector, sizeRows, sizeCols }
}

