# Security policy

`gueye` is a scientific Python package. The optional graphical interface starts an HTTP server bound to the loopback interface (`127.0.0.1`) and is intended for local interactive use only. It is not designed to be exposed to a network or used as a multi-user service.

Please report security concerns privately to the maintainer at `bira.gueye@ucad.edu.sn`, including the affected version and a minimal reproduction where possible. Do not include confidential data in public issues.

The interface applies request-size and grid-size limits. Users should still avoid opening the interface on untrusted machines or exposing its port through tunnels or proxies.
