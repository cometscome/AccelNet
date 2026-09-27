program runtime_probe
    use omp_lib
    implicit none
    integer :: x
    logical :: on_host
    on_host=.true.
    !$omp target map(from:on_host)
    on_host=omp_is_initial_device()
    !$omp end target
    if (on_host) error stop 'CPU fallback'
    x=1
    !$omp target enter data map(to:x)
    !$omp target map(alloc:x)
    x=x+1
    !$omp end target
    !$omp target update from(x)
    !$omp target exit data map(delete:x)
    if (x /= 2) error stop 'bad result'
    print *, 'runtime probe passed'
end program
