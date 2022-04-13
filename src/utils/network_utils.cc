#include "utils/network_utils.h"

#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/ip.h>
#include <sys/socket.h>
#include <sys/types.h>

#include "glog/logging.h"

namespace circinus {

int GetAvailablePort(int default_port) {
  struct sockaddr_in addr;
  addr.sin_port = htons(default_port);       // 0 means let system pick up an available port.
  addr.sin_family = AF_INET;                 // IPV4
  addr.sin_addr.s_addr = htonl(INADDR_ANY);  // set addr to any interface

  int sock = socket(AF_INET, SOCK_STREAM, 0);
  if (0 != bind(sock, (struct sockaddr*)&addr, sizeof(struct sockaddr_in))) {
    DLOG(WARNING) << "bind()";
    return 0;
  }
  socklen_t addr_len = sizeof(struct sockaddr_in);
  if (0 != getsockname(sock, (struct sockaddr*)&addr, &addr_len)) {
    DLOG(WARNING) << "getsockname()";
    return 0;
  }

  int ret = ntohs(addr.sin_port);
  close(sock);
  return ret;
}

int GetAvailablePort(int default_port, int n_tries) {
  auto target_port = default_port;
  auto port = 0;
  LOG(INFO) << "Try binding to " << target_port << "...";
  while ((port = GetAvailablePort(target_port)) == 0) {
    if (++target_port > default_port + n_tries) {
      LOG(FATAL) << "Cannot find an available port after " << n_tries << " probes.";
    }
    LOG(INFO) << "Can not bind to " << (target_port - 1) << ". Retry binding to " << target_port << "...";
  }
  return port;
}

}  // namespace circinus
