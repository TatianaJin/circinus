#pragma once

namespace circinus {

/** Get a free port on IPV4. */
int GetAvailablePort(int default_port = 0);

int GetAvailablePort(int default_port, int n_tries);

}  // namespace circinus
