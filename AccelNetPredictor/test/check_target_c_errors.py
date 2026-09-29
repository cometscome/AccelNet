"""Bad n2p2 files must return diagnostics without killing the embedding process."""
import ctypes as c
from pathlib import Path
import shutil
import sys
import tempfile

lib = c.CDLL(sys.argv[1])
data = Path(sys.argv[2])
create = lib.accelnet_target_create_n2p2
create.argtypes = [c.c_char_p, c.c_int, c.c_int, c.POINTER(c.c_int),
                   c.POINTER(c.c_double), c.POINTER(c.c_void_p), c.c_char_p]
create.restype = c.c_int
lib.accelnet_target_destroy.argtypes = [c.c_void_p]
lib.accelnet_target_get_species.argtypes = [c.c_void_p, c.c_int, c.c_int, c.c_char_p, c.c_char_p]


def load(path):
    count, cutoff, handle = c.c_int(99), c.c_double(99), c.c_void_p()
    error = c.create_string_buffer(512)
    status = create(str(path).encode(), -1, 0, c.byref(count), c.byref(cutoff), c.byref(handle), error)
    return status, count.value, cutoff.value, handle, error.value


base = (data / 'n2p2/input.nn').read_text()
angular = (data / 'n2p2-virial-angular/input.nn').read_text()
cases = {
    'empty': ('n2p2', '', None),
    'bad_count': ('n2p2', base.replace('number_of_elements 1', 'number_of_elements bad'), None),
    'huge_count': ('n2p2', base.replace('number_of_elements 1', 'number_of_elements 2147483647'), None),
    'duplicate_count': ('n2p2', base + '\nnumber_of_elements 1\n', None),
    'duplicate_elements': ('n2p2', base + '\nelements H\n', None),
    'missing_elements': ('n2p2', base.replace('elements H\n', ''), None),
    'unknown_element': ('n2p2', base.replace('elements H', 'elements Xx'), None),
    'duplicate_species': ('n2p2-virial-angular', angular.replace('elements H O', 'elements H H'), None),
    'type4g': ('n2p2', base.replace('2G-HDNNP', '4G-HDNNP'), None),
    'cutoff': ('n2p2', base.replace('cutoff_type 7 0.2', 'cutoff_type 99'), None),
    'alpha': ('n2p2', base.replace('cutoff_type 7 0.2', 'cutoff_type 9 0'), None),
    'normalization': ('n2p2', base.replace('conv_energy 2.0\n', ''), None),
    'nan_length': ('n2p2', base.replace('conv_length 3.0', 'conv_length NaN'), None),
    'layers': ('n2p2', base.replace('global_hidden_layers_short 0', 'global_hidden_layers_short 64'), None),
    'nodes': ('n2p2', base.replace('global_nodes_short', 'global_nodes_short bad'), None),
    'activation': ('n2p2', base.replace('global_activation_short l', 'global_activation_short x'), None),
    'activation_token': ('n2p2', base.replace('global_activation_short l', 'global_activation_short linear'), None),
    'radial': ('n2p2', base.replace('1.0 0.0 3.0', '1.0 0.0 -3.0'), None),
    'sf_type': ('n2p2', base.replace('symfunction_short H 2', 'symfunction_short H 99'), None),
    'sf_parameters': ('n2p2', base.replace('1.0 0.0 3.0', 'bad'), None),
    'unknown_neighbor': ('n2p2', base.replace('H 2 H', 'H 2 He'), None),
    'nan_radial': ('n2p2', base.replace('1.0 0.0 3.0', 'NaN 0.0 3.0'), None),
    'compact': ('n2p2', base.replace('H 2 H 1.0 0.0 3.0', 'H 20 H 0.0 3.0 bad'), None),
    'compact_interval': ('n2p2', base.replace('H 2 H 1.0 0.0 3.0', 'H 21 H H 0 3 30 10 p2'), None),
    'angular': ('n2p2', base.replace('H 2 H 1.0 0.0 3.0', 'H 3 H H 0.1 2 2 3.0'), None),
    'no_weights': ('n2p2', base, ('weights.001.data', None)),
    'short_weights': ('n2p2', base, ('weights.001.data', '0.1\n')),
    'bad_weights': ('n2p2', base, ('weights.001.data', 'bad\n')),
    'nan_weights': ('n2p2', base, ('weights.001.data', 'NaN\n0.1\n')),
    'no_scaling': ('n2p2-virial-angular', angular, ('scaling.data', None)),
    'bad_scaling': ('n2p2-virial-angular', angular, ('scaling.data', 'bad\n')),
    'bad_index': ('n2p2-virial-angular', angular, ('scaling.data', '1 -1 0 1 0 1\n')),
    'huge_index': ('n2p2-virial-angular', angular, ('scaling.data', '1 99999 0 1 0 1\n')),
    'nan_scaling': ('n2p2-virial-angular', angular, ('scaling.data', '1 1 0 NaN 0 1\n')),
    'zero_range': ('n2p2-virial-angular', angular, ('scaling.data', '1 1 0 0 0 1\n')),
    'scaling_interval': ('n2p2-virial-angular', angular.replace('scale_max_short 1.0', 'scale_max_short -1.0'), None),
}

status, count, cutoff, original, error = load(data / 'n2p2')
assert status == 0 and original.value and count == 1 and cutoff > 0, error
try:
    with tempfile.TemporaryDirectory(prefix='accelnet-c-errors-') as directory:
        root = Path(directory)
        for name, (fixture, text, edit) in cases.items():
            target = root / name
            shutil.copytree(data / fixture, target)
            (target / 'input.nn').write_text(text)
            if edit:
                file, value = edit
                if value is None:
                    (target / file).unlink()
                else:
                    (target / file).write_text(value)
            # Repeating failures also exercises file closure and partial-model cleanup.
            for _ in range(4):
                status, count, cutoff, failed, error = load(target)
                assert status != 0 and not failed.value and count == 0 and cutoff == 0 and error, (name, error)
            # A valid model can still be loaded after every failure.
            status, _, _, recovered, error = load(data / fixture)
            assert status == 0 and recovered.value, (name, error)
            lib.accelnet_target_destroy(recovered)
        symbol, error = c.create_string_buffer(17), c.create_string_buffer(512)
        assert lib.accelnet_target_get_species(original, 1, 17, symbol, error) == 0
        assert symbol.value == b'H'
finally:
    lib.accelnet_target_destroy(original)
print(f'{len(cases)} invalid n2p2 cases: diagnostics, repeated recovery and live-handle isolation passed')
